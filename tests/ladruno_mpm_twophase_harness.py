"""WP-232l harness: two-phase MPM particles against 232l-V's records (T0 i-ii, T1).

Reads tests/data/mpm/twophase/t0_oracle.json and t1_fe.json and the decks of
tests/mpm_twophase_cases.txt (blocks t0_modeA, t1_ref_fe_identity). It never reads
tests/mpm_twophase_reference.py (independence, wp-briefs "232l").

Run standalone to write a results record:
    python tests/ladruno_mpm_twophase_harness.py --bin <dist\\bin> [--out <json>]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "mpm" / "twophase"
T0_REC = DATA / "t0_oracle.json"
T1_REC = DATA / "t1_fe.json"

T0_TOL = 1e-9          # band: max |C++ - oracle| / max |oracle|, per field
T1_A_TOL = 1e-12       # Frobenius, relative
T1_B_TOL = 1e-12       # 2-norm, relative
T1_STEP_TOL = 1e-10    # first step u, p, relative


def as_list(v):
    if isinstance(v, (list, tuple)):
        return [float(x) for x in v]
    return [float(v)]


def rel_max(a, b):
    """max |a - b| / max |b| over flattened arrays (0 when both vanish)."""
    fa, fb = _flat(a), _flat(b)
    if len(fa) != len(fb):
        return float("inf")
    den = max((abs(x) for x in fb), default=0.0)
    num = max((abs(x - y) for x, y in zip(fa, fb)), default=0.0)
    if den == 0.0:
        return num
    return num / den


def _flat(x):
    if isinstance(x, (list, tuple)):
        out = []
        for y in x:
            out.extend(_flat(y))
        return out
    return [float(x)]


def _key(c, nd):
    return tuple(round(float(c[i]), 9) for i in range(nd))


# ---------------------------------------------------------------------------
#  T0: mode A, route E (case t0_modeA)
# ---------------------------------------------------------------------------
def build_t0_modeA(ops, rec):
    c = rec["case_modeA"]
    fl = c["fluid"]
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, c["E"], c["nu"], c["rho"])
    ops.mpm("grid", 232, "-origin", c["x0"][0], c["x0"][1], "-h", c["h"], "-cells", c["cells"][0],
            c["cells"][1], "-basis", "gimp")
    ops.mpm("mode", "implicit")
    ops.mpm("kinematics", "small")
    ops.mpm("body", 1, "box", c["box"][0][0], c["box"][0][1], c["box"][1][0], c["box"][1][1],
            "-ppc", c["ppc"], "-mat", 1, "-rho", c["rho"])
    stab = fl["stab"]
    stab_args = ["-stab", "off"] if stab[0] == "off" else ["-stab", "auto", float(stab[1])]
    ops.mpm("twoPhase", 1, "-Kf", fl["Kf"], "-poro", fl["poro"], "-rhoF", fl["rhoF"],
            "-perm", fl["perm"][0], fl["perm"][1], *stab_args)
    ops.mpm("gravity", c["gravity"][0], c["gravity"][1])
    for face, comps in c["fix"]:
        ops.mpm("fix", face, "".join("xyz"[i] for i in comps))
    for face, p0 in c["pfix"]:
        ops.mpm("pfix", face, float(p0))
    ops.mpm("rayleigh", c["rayleigh"][0], c["rayleigh"][1])
    ops.mpm("flip", c["flip"])
    ops.mpm("pressureBlend", c["pblend"])
    h0 = rec["t0_ii"]["history"][0]["particles"]
    ops.mpm("particleVelocity", *_flat(h0["v"]))
    ops.mpm("particleAcceleration", *_flat(h0["a"]))
    ops.mpm("particlePressure", *(list(h0["p"]) + list(h0["pd"]) + list(h0["pdd"])))
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.algorithm("Linear")
    ops.integrator("Newmark", c["newmark"][0], c["newmark"][1])
    ops.analysis("Transient")


def particle_fields(ops, nd=2, nv=3):
    def grab(f, w):
        v = as_list(ops.mpm("particleData", f))
        return [v[i:i + w] for i in range(0, len(v), w)] if w > 1 else v
    return {"x": grab("x", nd), "v": grab("v", nd), "a": grab("a", nd), "sig": grab("sig", nv),
            "p": grab("p", 1), "pd": grab("pd", 1), "pdd": grab("pdd", 1)}


def run_t0_i(ops, rec, tmp):
    """T0 (i): each listed particle's blocks against the oracle's (§7a Q1, Q2)."""
    s = rec["case_modeA"]["dump"]["s"]
    out = []
    for e in rec["t0_i"]:
        build_t0_modeA(ops, rec)
        path = str(Path(tmp) / ("blocks_%d.json" % e["particle"]))
        ops.mpm("debug", "particleBlocks", e["particle"], path, s["u"], s["p"], s["v"], s["pd"],
                s["a"], s["pdd"])
        d = json.loads(Path(path).read_text())
        nd, ndf = d["ndm"], d["ndf"]
        mine = {_key(c, nd): k for k, c in enumerate(d["coords"])}
        perm = [mine.get(_key(c, nd), -1) for c in e["coords"]]
        row = {"particle": e["particle"], "nen": d["nen"], "nen_oracle": len(e["coords"])}
        if -1 in perm or len(perm) != d["nen"]:
            row["error"] = "support nodes differ from the oracle's"
            out.append(row)
            continue
        dof = [ndf * perm[k] + c for k in range(len(perm)) for c in range(ndf)]
        for f in ("K", "C", "M"):
            A = d[f]
            row[f] = rel_max([[A[i][j] for j in dof] for i in dof], e[f])
        for f in ("R", "R_static"):
            row[f] = rel_max([d[f][i] for i in dof], e[f])
        row["S"] = rel_max([d["S"][k] for k in perm], e["S"])
        row["dS"] = rel_max([d["dS"][k] for k in perm], e["dS"])
        row["alpha_st"] = rel_max([d["alpha_st"]], [e["alpha_st"]])
        row["max"] = max(v for k, v in row.items() if k in ("K", "C", "M", "R", "R_static", "S", "dS",
                                                             "alpha_st"))
        out.append(row)
    return out


def run_t0_ii(ops, rec):
    """T0 (ii): three mode-A Newmark steps; every particle field after each (history[k])."""
    build_t0_modeA(ops, rec)
    hist = rec["t0_ii"]["history"]
    dt = rec["case_modeA"]["dt"]
    rows = []
    for k, h in enumerate(hist):
        if k > 0:
            ops.mpm("analyze", 1, dt)
        mine = particle_fields(ops)
        row = {"step": k}
        for f in ("x", "v", "a", "sig", "p", "pd", "pdd"):
            if f in h["particles"]:
                row[f] = rel_max(mine[f], h["particles"][f])
        row["max"] = max(v for kk, v in row.items() if kk != "step")
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
#  T1: FE identity (case t1_ref_fe_identity)
# ---------------------------------------------------------------------------
def build_t1(ops, rec, leg):
    L = rec["legs"][leg]
    c = rec["case"]
    nd = L["ndm"]
    mesh = L["mesh"]
    ncell = {"q4_1": [1, 1], "q4_2x2": [2, 2], "h8_1": [1, 1, 1]}[mesh]
    h = c["h"]
    ops.wipe()
    ops.model("basic", "-ndm", nd, "-ndf", nd)
    if L["material"] == "elastic":
        ops.nDMaterial("ElasticIsotropic", 1, c["E"], c["nu"], c["rho"])
    else:
        ops.nDMaterial("DruckerPrager", 1, *[float(x) for x in c["dp_args"]])
    ops.mpm("grid", 232, "-origin", *([0.0] * nd), "-h", h, "-cells", *ncell, "-basis", "linear")
    ops.mpm("mode", "implicit")
    ops.mpm("kinematics", "small")
    lp = [h / 4.0] * nd
    for i, p in enumerate(L["particles"]):
        rho = p["m"] / p["V"]
        ops.mpm("particle", i, 1, 1, rho, p["V"], "-x", *p["x"], "-lp", *lp)
    stab = L["stab"]                         # "off" or ["auto", a0]
    stab_args = ["-stab", "off"] if stab == "off" else ["-stab", "auto", float(stab[1])]
    ops.mpm("twoPhase", 1, "-Kf", c["Kf"], "-poro", c["poro"], "-rhoF", c["rhoF"],
            "-perm", *([c["kbar"]] * nd), *stab_args)
    g = [0.0] * nd
    g[1] = -c["g"]
    ops.mpm("gravity", *g)
    P = L["particles"]
    ops.mpm("particleVelocity", *[x for p in P for x in p["v"]])
    ops.mpm("particleAcceleration", *[x for p in P for x in p["a"]])
    ops.mpm("particlePressure", *([p["p"] for p in P] + [p["pd"] for p in P] + [p["pdd"] for p in P]))
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.algorithm("Linear")
    ops.integrator("Newmark", c["gamma"], c["beta"])
    ops.analysis("Transient")
    return L


def run_t1(ops, rec, tmp, legs=None):
    out = []
    dt = rec["case"]["dt"]
    for leg in (legs or list(rec["legs"].keys())):
        L = build_t1(ops, rec, leg)
        nd = L["ndm"]
        path = str(Path(tmp) / ("assemble_%s.json" % leg))
        ops.mpm("debug", "assemble", path, dt)
        d = json.loads(Path(path).read_text())
        row = {"leg": leg, "n": d["n"], "n_fe": len(L["dofs"])}
        mine = {(_key(lab[1:1 + nd], nd), int(lab[1 + nd])): e for e, lab in enumerate(d["dofs"])}
        idx = [mine.get((_key(lab[1:1 + nd], nd), int(lab[1 + nd])), -1) for lab in L["dofs"]]
        if -1 in idx or d["n"] != len(L["dofs"]):
            row["error"] = "equation sets differ"
            out.append(row)
            continue
        A = [[d["A"][i][j] for j in idx] for i in idx]
        b = [d["b"][i] for i in idx]
        num = math.sqrt(sum((A[i][j] - L["A"][i][j]) ** 2 for i in range(len(idx)) for j in range(len(idx))))
        den = math.sqrt(sum(x * x for r in L["A"] for x in r))
        row["A_rel"] = num / den
        num = math.sqrt(sum((x - y) ** 2 for x, y in zip(b, L["b"])))
        den = math.sqrt(sum(y * y for y in L["b"]))
        row["b_rel"] = num / den if den > 0 else num
        # the first step, from a fresh build (assemble left the step unrun)
        build_t1(ops, rec, leg)
        ops.mpm("analyze", 1, dt)
        tags = [int(t) for t in as_list(ops.mpm("info", "gridNodeTags"))]
        bypos = {}
        for t in tags:
            bypos[_key(ops.nodeCoord(t), nd)] = t
        u_m, u_f, p_m, p_f = [], [], [], []
        miss = 0
        for nrec in L["nodes"]:
            t = bypos.get(_key(nrec["x"], nd))
            if t is None:
                miss += 1
                continue
            dsp = as_list(ops.nodeDisp(t))
            u_m.extend(dsp[:nd])
            u_f.extend(nrec["u1"])
            p_m.append(dsp[nd])
            p_f.append(nrec["p1"])
        row["missing_nodes"] = miss
        row["u_rel"] = rel_max(u_m, u_f)
        row["p_rel"] = rel_max(p_m, p_f)
        out.append(row)
    return out


def t1_pass(row):
    return ("error" not in row and row["missing_nodes"] == 0 and row["A_rel"] <= T1_A_TOL and
            row["b_rel"] <= T1_B_TOL and row["u_rel"] <= T1_STEP_TOL and row["p_rel"] <= T1_STEP_TOL)


def load(path):
    return json.loads(Path(path).read_text())


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="dist\\bin holding opensees.pyd")
    ap.add_argument("--out", default=str(DATA / "t0_t1_mpm.json"))
    a = ap.parse_args(argv)
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(a.bin)
    sys.path.insert(0, a.bin)
    import opensees as ops  # noqa: E402
    build = ops.ladrunoBuild() if hasattr(ops, "ladrunoBuild") else "unknown"
    r0, r1 = load(T0_REC), load(T1_REC)
    with tempfile.TemporaryDirectory() as tmp:
        res = {"leg": "t0_t1_mpm", "build": build, "bands": {"t0": T0_TOL, "t1_A": T1_A_TOL,
               "t1_b": T1_B_TOL, "t1_step": T1_STEP_TOL},
               "t0_i": run_t0_i(ops, r0, tmp), "t0_ii": run_t0_ii(ops, r0), "t1": run_t1(ops, r1, tmp)}
    res["t0_i_pass"] = all("error" not in r and r["max"] <= T0_TOL for r in res["t0_i"])
    res["t0_ii_pass"] = all(r["max"] <= T0_TOL for r in res["t0_ii"])
    res["t1_pass"] = all(t1_pass(r) for r in res["t1"])
    Path(a.out).write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ("build", "t0_i_pass", "t0_ii_pass", "t1_pass")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
