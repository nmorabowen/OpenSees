"""WP-232k harness, rows S1-S6 and S10 (interface tests/mpm_soil_interface.md; gates in
wp-briefs.md "232k" and its freeze). Written from the interface and the reference blocks of
tests/mpm_soil_cases.txt only; never reads tests/mpm_soil_reference.py. Units kN, t, m, s, kPa.

S1  one-point identity: a 2 x 2 (x 2)-cell block, ppc 2, every grid node with mass prescribed
    the affine velocity v = L (x - ref) (`mpm gridvel ... affine`, `setAffine`), L = (dF - I)/dt
    from the record's per-step dF, so every particle sees dF exactly (linear completeness;
    `mpm massCutoff 0` keeps every node with mass). Mode A: zero free equations, as the FE
    driver. Mode F: dt 1e-4. Compared per step with GP 1 of the record: the material-frame
    stress (`particleData sigHat`; 3D Voigt; plane strain sigma_zz from `particleData sigZZ`) and
    the material state vector.
    Tolerance (freeze): rel 1e-8 before fe_self_consistency.first_step_gp_rel_spread_gt_1e8,
    rel 1e-3 from it on (never: 1e-8 throughout).
S2  the K0 column (fixed base, `fix ymin y`) + a surcharge rho g l0 on the top particle row, as a
    per-particle body force ramped over 10 steps after the 2-step hold (mode A).
S4  the low-g twin (g 0.5, l0 2): the flip map per particle in both modes against the record's
    rule (p < Pmin -> sigma = -Pmin I, alpha = 0; else sigma bit for bit).
S10 classic Tcl against Python on one S1 leg.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "tests" / "data" / "mpm" / "soil"
HOST = 232
TS_TAG = 9001
RHO, G = 2.0, 9.81


def _l(v):
    return [float(x) for x in (v if isinstance(v, (list, tuple)) else [v])]


def _num(x):
    try:
        return int(x)
    except ValueError:
        try:
            return float(x)
        except ValueError:
            return x


def deck_args(cmd, tag=1, extra=()):
    """'nDMaterial LadrunoHS 1 a b ...' -> ['LadrunoHS', tag, a, b, ...]."""
    t = cmd.split()
    assert t[0] == "nDMaterial"
    return [t[1], tag] + [_num(x) for x in t[3:]] + list(extra)


def load(name):
    return json.loads((DATA / name).read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- #
#  S1                                                                          #
# --------------------------------------------------------------------------- #
def s1_tolerance(leg_rec, step):
    """The freeze's rule: 1e-8 before the FE's own Gauss-point spread exceeds 1e-8, 1e-3 after."""
    k = (leg_rec.get("fe_self_consistency") or {}).get("first_step_gp_rel_spread_gt_1e8")
    return 1e-8 if (k is None or step < k) else 1e-3


def s1_commands(rec, leg, mode, dt):
    """[(command, args)] that build the S1 block (no steps)."""
    L = rec["legs"][leg]
    nd = 3 if L["dim"] == "3d" else 2
    mat = L["material"]
    cmds = [("wipe", []), ("model", ["basic", "-ndm", nd, "-ndf", nd]),
            ("nDMaterial", deck_args(rec["decks"][mat]["cmd"])),
            ("mpm", ["grid", HOST, "-origin", *([-1.0] * nd), "-h", 0.5, "-cells", *([6] * nd)]),
            ("mpm", ["mode", "quasiStatic" if mode == "A" else "explicit"]),
            ("mpm", ["body", 1, "box", *([0.0] * nd), *([1.0] * nd), "-ppc", 2, "-mat", 1,
                     "-rho", RHO]),
            ("mpm", ["massCutoff", 0.0]),
            ("mpm", ["gridvel", "all", "box", *([-1.6] * nd), *([2.6] * nd), "affine",
                     *([0.0] * (nd * nd)), "ref", *([0.5] * nd)])]
    if mode == "A":
        cmds += [("constraints", ["Transformation"]), ("numberer", ["RCM"]),
                 ("system", ["UmfPack"]), ("test", ["NormDispIncr", 1e-12, 10]),
                 ("algorithm", ["Newton"]), ("integrator", ["LoadControl", dt]),
                 ("analysis", ["Static"])]
    cmds.append(("mpm", ["stage", 1, 0]))
    return cmds, nd


def s1_L(rec, leg, k, nd, dt):
    dF = rec["legs"][leg]["dF"][k]
    return [(dF[a][b] - (1.0 if a == b else 0.0)) / dt for a in range(nd) for b in range(nd)]


def _sig6(sig, szz, nd, p):
    if nd == 3:
        return sig[6 * p:6 * p + 6]
    s = sig[3 * p:3 * p + 3]
    return [s[0], s[1], szz[p], s[2], 0.0, 0.0]


def s1_run(ops, leg, mode, rec=None, nsteps=None):
    rec = rec or load("s1_onepoint.json")
    L = rec["legs"][leg]
    dt = 1.0 if mode == "A" else 1e-4
    cmds, nd = s1_commands(rec, leg, mode, dt)
    for name, args in cmds:
        getattr(ops, name)(*args)
    n = len(L["dF"]) if nsteps is None else nsteps
    flip = int(L["flip_before_step"])
    out = {"leg": leg, "mode": mode, "steps": [], "first_substep_diff": None}
    for k in range(n):
        if k + 1 == flip:
            ops.mpm("stage", 1, 1)
        ops.mpm("gridvel", "all", "setAffine", *s1_L(rec, leg, k, nd, dt))
        ops.mpm("analyze", 1, dt)
        # the record's sig is GP 1's MATERIAL-frame stress (the unrotated sigma_hat; on the shear
        # legs it differs from the rotated Cauchy stress by the rotation, ~1e-3): compare sigHat
        sig = _l(ops.mpm("particleData", "sigHat"))
        szz = _l(ops.mpm("particleData", "sigZZ"))
        st = _l(ops.mpm("particleResponse", "state"))
        N = len(szz)
        ns = len(st) // N
        ref = L["steps"][k]
        rs, rst = ref["sig"], ref["state"]
        scale = max(abs(v) for v in rs) or 1.0
        sscale = max(abs(v) for v in rst) or 1.0
        es = max(max(abs(a - b) for a, b in zip(_sig6(sig, szz, nd, p), rs)) for p in range(N)) / scale
        et = max(max(abs(st[ns * p + i] - rst[i]) for i in range(min(ns, len(rst))))
                 for p in range(N)) / sscale
        row = {"step": k + 1, "err_sig": es, "err_state": et, "tol": s1_tolerance(L, k + 1)}
        if L["material"] == "sanisand" and "substeps" in ref:
            sub = _l(ops.mpm("particleResponse", "substeps"))
            mine = sub[0::len(sub) // N] if N else []
            row["substeps"] = mine[0] if mine else None
            if out["first_substep_diff"] is None and mine and mine[0] != ref["substeps"][0]:
                out["first_substep_diff"] = k + 1
        out["steps"].append(row)
    out["max_err"] = max(max(r["err_sig"], r["err_state"]) for r in out["steps"])
    out["passed"] = all(r["err_sig"] <= r["tol"] and r["err_state"] <= r["tol"] for r in out["steps"])
    return out


def s1_tcl_deck(rec, leg, nsteps=None):
    """The same S1 leg as a classic-Tcl deck (mode A); prints the final sigma list."""
    L = rec["legs"][leg]
    cmds, nd = s1_commands(rec, leg, "A", 1.0)
    lines = []
    for name, args in cmds:
        lines.append(" ".join([name] + [repr(a) if isinstance(a, float) else str(a) for a in args]))
    n = len(L["dF"]) if nsteps is None else nsteps
    flip = int(L["flip_before_step"])
    for k in range(n):
        if k + 1 == flip:
            lines.append("mpm stage 1 1")
        Lk = " ".join(repr(v) for v in s1_L(rec, leg, k, nd, 1.0))
        lines.append(f"mpm gridvel all setAffine {Lk}")
        lines.append("mpm analyze 1 1.0")
    lines.append('puts "SIG [mpm particleData sig]"')
    return "\n".join(lines) + "\n"


def run_tcl(deck, exe, timeout=600):
    with tempfile.TemporaryDirectory() as td:
        f = Path(td) / "s1.tcl"
        f.write_text(deck, encoding="utf-8")
        env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
        r = subprocess.run([exe, str(f)], capture_output=True, text=True, cwd=td,
                           timeout=timeout, env=env)
    for line in (r.stdout + r.stderr).splitlines():
        if line.startswith("SIG "):
            return [float(x) for x in line[4:].split()]
    raise RuntimeError("no SIG line from the Tcl run:\n" + r.stdout[-2000:] + r.stderr[-2000:])


# --------------------------------------------------------------------------- #
#  S2 / S4 -- the column with a fixed base                                     #
# --------------------------------------------------------------------------- #
CAP_E, CAP_NU, CAP_RHO = 214234.27337502, 0.3129, RHO
CAP_T = 0.5                                      # the cap thickness (S2a)
CAP_SOIL = False                                 # S2a variant: the cap is the soil material itself


def column_build(ops, mat_cmd, basis, mode, l0=10.0, g=G, s0=0, cap=False):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial(*deck_args(mat_cmd))
    cells_y = int(round(l0 / 0.5)) + 1 + (int(round(CAP_T / 0.5)) if cap else 0)
    ops.mpm("grid", HOST, "-origin", 0.0, 0.0, "-h", 0.5, "-cells", 2, cells_y, "-basis", basis)
    ops.mpm("mode", "quasiStatic" if mode == "A" else "explicit")
    ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, l0, "-ppc", 2, "-mat", 1, "-rho", RHO)
    if cap:
        # S2's surcharge carrier (amendment S2a, 232k-C): a bonded elastic cap one cell thick,
        # the soil's density (so its nodes stay active under the mode-A mass cutoff); its weight
        # is cancelled by its body force, and the surcharge is that body force, so it reaches
        # the soil at its top surface y = l0, as the FE's nodal top loads do. Not measured.
        if not CAP_SOIL:
            ops.nDMaterial("ElasticIsotropic", 98, CAP_E, CAP_NU)
        ops.mpm("body", 3, "box", 0.0, l0, 1.0, l0 + CAP_T, "-ppc", 2, "-mat",
                1 if CAP_SOIL else 98, "-rho", CAP_RHO)
    ops.mpm("fix", "xmin", "x")
    ops.mpm("fix", "xmax", "x")
    ops.mpm("fix", "ymin", "y")
    ops.mpm("gravity", 0.0, -g)
    ops.timeSeries("Path", TS_TAG, "-time", 0.0, 1.0, 1.0e9, "-values", 0.0, 1.0, 1.0)
    ops.mpm("timeSeries", TS_TAG)
    if mode == "A":
        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1.0e-10, 25)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 0.1)
        ops.analysis("Static")
    else:
        ops.mpm("damping", "local", 0.8)
    ops.mpm("stage", 1, s0)


def _column_state(ops):
    x = _l(ops.mpm("particleData", "x"))
    return dict(x=x, sig=_l(ops.mpm("particleData", "sig")),
                szz=_l(ops.mpm("particleData", "sigZZ")))


def _column_errors(st, Y, rows, rho_g_l0):
    """max |syy - ref| and RMS (sxx - ref) over particles, / rho g l0."""
    N = len(Y)
    e_yy, sq = 0.0, 0.0
    for k in range(N):
        r = min(rows, key=lambda q: abs(q["Y"] - Y[k]))
        e_yy = max(e_yy, abs(st["sig"][3 * k + 1] - r["syy"]))
        sq += (st["sig"][3 * k] - r["sxx"]) ** 2
    return e_yy / rho_g_l0, math.sqrt(sq / N) / rho_g_l0


def s2_run(ops, mat, basis, nsur=10, cap=True):
    """cap=True (amendment S2a): the surcharge rides on a nearly weightless bonded elastic cap
    above the column. cap=False: the 232f footing mechanism, a body force on the soil's top
    particle row -- the loaded row then carries about half the load in its own stress (the load
    acts at its centre, not at the surface), which the max-over-particles syy gate reads as an
    error of ~q/2 (diagnosis record)."""
    rec = load("s2_column.json")
    leg = rec["legs"][f"{mat}_nsur{nsur}"]
    column_build(ops, rec["decks"][mat]["cmd"], basis, "A", cap=cap)
    x0 = _l(ops.mpm("particleData", "x"))
    body = _l(ops.mpm("particleData", "body"))
    Nall = len(body)
    N = sum(1 for b in body if b == 1.0)
    Y = [x0[2 * k + 1] for k in range(N)]
    rgl = RHO * G * 10.0
    out = {"mat": mat, "basis": basis, "checkpoints": {}}

    def check(name):
        st = _column_state(ops)
        st = {"sig": st["sig"][:3 * N], "szz": st["szz"][:N], "x": st["x"][:2 * N]}
        ey, ex = _column_errors(st, Y, leg[name]["particles"], rgl)
        out["checkpoints"][name] = {"syy_max": ey, "sxx_rms": ex}

    m = _l(ops.mpm("particleData", "m"))
    ytop = max(Y)
    q = float(leg.get("surcharge", rgl))
    carriers = [k for k in range(Nall) if body[k] == 3.0] if cap else \
        [k for k in range(N) if abs(Y[k] - ytop) < 1e-9]

    def step(lam, frac):
        """one quasi-static step; the cap's weight cancelled at the step-end load factor lam,
        and the surcharge fraction frac on the carriers (per unit mass)"""
        bf = []
        for k in range(Nall):
            b = 0.0
            if cap and body[k] == 3.0:
                b += lam * G
            if k in carriers:
                b += -(q / len(carriers)) * frac / m[k]
            bf += [0.0, b]
        ops.mpm("particleBodyForce", *bf)
        ops.mpm("analyze", 1, 0.1)

    for i in range(1, 11):
        step(i / 10.0, 0.0)
    check("end_stage0")
    ops.mpm("stage", 1, 1)
    check("after_flip")
    step(1.0, 0.0)
    step(1.0, 0.0)
    check("after_hold")
    # the surcharge q = rho g l0 (per unit width W = 1), ramped over nsur steps
    for i in range(1, nsur + 1):
        step(1.0, i / nsur)
    check("after_surcharge")
    out["cap"] = cap
    out["passed"] = all(c["syy_max"] <= 0.03 and c["sxx_rms"] <= 0.02
                        for n, c in out["checkpoints"].items() if n != "after_hold")
    return out


def s4_flip_map(ops, mode, basis="gimp", g=0.5, l0=2.0):
    """Stage 0 under low g, then the walk; returns per-particle (pre, post) and the check of the
    record's rule: p < Pmin -> sigma = -Pmin I and alpha = 0, else sigma bit for bit; alpha_in =
    alpha after the flip (-flipAlphaIn init)."""
    rec = load("s4_lowg.json")
    pmin = float(rec["Pmin"])
    column_build(ops, rec["decks"]["sanisand"]["cmd"], basis, mode, l0=l0, g=g)
    if mode == "A":
        ops.mpm("analyze", 10, 0.1)
    else:
        ops.mpm("analyze", 3000, 5e-4)
    pre = _column_state(ops)
    pst = _l(ops.mpm("particleResponse", "state"))
    ops.mpm("stage", 1, 1)
    post = _column_state(ops)
    qst = _l(ops.mpm("particleResponse", "state"))
    N = len(pre["szz"])
    reset, bad = 0, []
    for k in range(N):
        s = pre["sig"][3 * k:3 * k + 3]
        p = -(s[0] + s[1] + pre["szz"][k]) / 3.0
        t = post["sig"][3 * k:3 * k + 3] + [post["szz"][k]]
        a = qst[26 * k + 6:26 * k + 12]
        ain = qst[26 * k + 18:26 * k + 24]
        if p < pmin:
            reset += 1
            ok = t == [-pmin, -pmin, 0.0, -pmin] or (t[0] == -pmin and t[1] == -pmin and
                                                    t[2] == 0.0 and t[3] == -pmin)
            ok = ok and all(v == 0.0 for v in a)
        else:
            ok = t[:3] == s and t[3] == pre["szz"][k]
        ok = ok and ain == a
        if not ok:
            bad.append(k)
    return {"mode": mode, "particles": N, "reset": reset, "bad": bad,
            "passed": reset >= 1 and not bad, "pre_state": pst[:26], "post_state": qst[:26]}


# --------------------------------------------------------------------------- #
#  S4, the same-state half (orchestrator's ruling at the S4 freeze)            #
# --------------------------------------------------------------------------- #
def s4_affine(ops, mode, row, n=10):
    """The S1 block (SANISAND, plane strain), an affine oedometric stage-0 path to the strain of
    the S4 record's row `row` (its pre-flip elastic strain eps_yy, state[1], compression
    positive), then `mpm stage 1 1`. Returns the particles' pre- and post-flip sigma (3D Voigt,
    sigma_zz from sigZZ) and SANISAND `state`, and the record row. Each of the n steps takes the
    same strain increment d = eps / n through the Hughes-Winget stretch lam = (2 - d)/(2 + d)
    (2 (lam - 1)/(lam + 1) = -d), so the hypo total strain is eps up to round-off."""
    rec = load("s4_lowg.json")
    r = rec["rows"][row]
    s1 = load("s1_onepoint.json")
    dt = 1.0 if mode == "A" else 1e-4
    cmds, nd = s1_commands(s1, "sanisand_ps_oedo", mode, dt)
    for name, args in cmds:
        getattr(ops, name)(*args)
    eps = r["pre_flip"]["state"][1]
    d = eps / n
    lam = (2.0 - d) / (2.0 + d)
    for _ in range(n):
        ops.mpm("gridvel", "all", "setAffine", 0.0, 0.0, 0.0, (lam - 1.0) / dt)
        ops.mpm("analyze", 1, dt)

    def snap():
        sig = _l(ops.mpm("particleData", "sig"))
        szz = _l(ops.mpm("particleData", "sigZZ"))
        N = len(szz)
        return {"sig": [_sig6(sig, szz, 2, p) for p in range(N)],
                "state": _l(ops.mpm("particleResponse", "state"))}

    pre = snap()
    ops.mpm("stage", 1, 1)
    post = snap()
    return {"mode": mode, "row": row, "below_Pmin": r["below_Pmin"], "pre": pre, "post": post,
            "ref_pre": r["pre_flip"], "ref_post": r["after_flip"]}


def rel_diff(a, b):
    sc = max(max(abs(v) for v in b), 1e-300)
    return max(abs(x - y) for x, y in zip(a, b)) / sc


# --------------------------------------------------------------------------- #
#  S3 -- the strip footing (232j L4 domain), mode A QS, K9's ladder            #
# --------------------------------------------------------------------------- #
LADDER = [("KrylovNewton", 1e-6, 40), ("NewtonLineSearch", 1e-6, 60), ("KrylovNewton", 1e-5, 80)]
RELAXED = 3            # the third try is the relaxed one
H_MIN = 5e-4           # the halving floor, of the 2 mm (pseudo dt 1) step
Q0 = 20.0              # the surcharge, kPa
B_HALF = 1.0


def s3_build(ops, mat_cmd, r, antilocking="fbar", extra_opts=()):
    h = 0.5 / 2 ** (r - 1)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial(*deck_args(mat_cmd, extra=extra_opts))
    ops.mpm("grid", HOST, "-origin", 0.0, -6.0, "-h", h, "-cells", int(round(12 / h)),
            int(round(7 / h)), "-basis", "gimp")
    ops.mpm("mode", "quasiStatic")
    ops.mpm("body", 1, "box", 0.0, -6.0, 12.0, 0.0, "-ppc", 2, "-mat", 1, "-rho", RHO)
    ops.mpm("pforce", 1, "box", 0.0, -h / 2.0, 12.0, 0.0, 0.0, -Q0 * h / 2.0)
    ops.mpm("fix", "xmin", "x")
    ops.mpm("fix", "xmax", "x")
    ops.mpm("fix", "ymin", "xy")
    ops.mpm("gridvel", "footing", "box", -h, 0.0, B_HALF, 1.0e9, "v", 0.0, -0.002,
            "start", 12.0, "track", "y")
    ops.mpm("gravity", 0.0, -G)
    ops.timeSeries("Path", TS_TAG, "-time", 0.0, 10.0, 1.0e9, "-values", 0.0, 1.0, 1.0)
    ops.mpm("timeSeries", TS_TAG)
    ops.mpm("antiLocking", antilocking)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", 1e-6, 40)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ops.mpm("cut", 2.0, 0)
    ops.mpm("stage", 1, 0)


def _ladder_step(ops, t_target, state):
    """One accepted step to t_target with K9's ladder; halving below the step, floor H_MIN."""
    h = state["h"]
    while True:
        t = float(ops.mpm("info", "time"))
        if abs(t - t_target) <= 1e-9 * max(1.0, abs(t_target)):
            state["h"] = 1.0
            return True
        hh = min(h, t_target - t)
        ok = None
        for i, (alg, tol, it) in enumerate(LADDER):
            ops.algorithm(alg)
            ops.test("NormUnbalance", tol, it)
            try:
                ops.mpm("analyze", 1, hh)
                ok = i
                break
            except Exception:
                state["failed"] += 1
        if ok is None:
            h = hh / 2.0
            state["halvings"] += 1
            if h < H_MIN:
                return False
            continue
        if ok + 1 >= RELAXED:
            state["relaxed"] += 1
        h = min(2.0 * hh, 1.0)


def s3_run(ops, mat_cmd, r, antilocking="fbar", extra_opts=(), s_end=0.1, log=None):
    import time as _t
    w0 = _t.time()
    s3_build(ops, mat_cmd, r, antilocking, extra_opts)
    st = {"h": 1.0, "failed": 0, "relaxed": 0, "halvings": 0}
    out = {"r": r, "antilocking": antilocking, "curve": [], "reached": False, "stopped": None}
    for k in range(1, 13):                          # stage 0 (10), the flip, 2 holds
        if k == 11:
            ops.mpm("stage", 1, 1)
        if not _ladder_step(ops, float(k), st):
            out["stopped"] = f"t = {k - 1} (staging)"
            break
    n_push = int(round(s_end / 0.002))
    if out["stopped"] is None:
        for k in range(1, n_push + 1):
            if not _ladder_step(ops, 12.0 + k, st):
                out["stopped"] = f"s = {0.002 * (k - 1):.4f}"
                break
            q = -_l(ops.mpm("condForce", "footing"))[1] / B_HALF
            out["curve"].append([0.002 * k, q])
            if log and k % 10 == 0:
                log(f"S3 r{r}: s = {0.002 * k:.3f}, q = {q:.2f}, failed {st['failed']}")
        else:
            out["reached"] = True
    out["q_at"] = {}
    for sb in (0.01, 0.02, 0.05):
        for si, q in out["curve"]:
            if abs(si - sb * 2.0) < 1e-9:
                out["q_at"][str(sb)] = q
    out.update(failed_attempts=st["failed"], relaxed=st["relaxed"], halvings=st["halvings"],
               wall_s=_t.time() - w0)
    return out


def s3_judge(out, ref):
    rows = {}
    for sb in ("0.01", "0.02", "0.05"):
        q, qr = out["q_at"].get(sb), ref["q_at"][sb]
        rows[sb] = {"mpm": q, "fe": qr, "ratio": (q / qr) if q is not None else None,
                    "pass": q is not None and abs(q / qr - 1.0) <= 0.10}
    return {"rows": rows, "reached": out["reached"],
            "passed": out["reached"] and all(v["pass"] for v in rows.values())}


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["s1", "s2", "s3"])
    ap.add_argument("--leg")
    ap.add_argument("--mode", default="A")
    ap.add_argument("--mat", default="hs")
    ap.add_argument("--basis", default="gimp")
    ap.add_argument("--r", type=int, default=1)
    ap.add_argument("--al", default="fbar")
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from _testbed import ops   # noqa: E402
    if a.what == "s1":
        res = s1_run(ops, a.leg, a.mode)
        name = f"s1_mpm_{a.leg}_{a.mode}.json"
    elif a.what == "s2":
        res = s2_run(ops, a.mat, a.basis)
        name = f"s2_mpm_{a.mat}_{a.basis}.json"
    else:
        ref = load(f"s3_fe_{a.mat}_n{4 * 2 ** (a.r - 1)}.json")
        res = s3_run(ops, ref["decks"][a.mat]["cmd"], a.r, a.al, log=print)
        res["judge"] = s3_judge(res, ref)
        name = f"s3_mpm_{a.mat}_r{a.r}_{a.al}.json"
    res["build"] = str(ops.ladrunoBuild())
    path = Path(a.out) if a.out else DATA / name
    path.write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in res.items() if k not in ("steps", "curve")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
