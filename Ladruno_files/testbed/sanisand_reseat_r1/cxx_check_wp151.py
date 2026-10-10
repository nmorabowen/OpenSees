import json, os, sys, math, statistics, collections
W = sys.argv[1]; BIN = os.path.join(W, "dist", "bin")
os.add_dll_directory(BIN); sys.path.insert(0, BIN); sys.path.insert(0, os.path.join(W, "tests"))
import opensees as ops
assert os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) == os.path.normcase(BIN), ops.__file__
import wp151_reseat_tools as T
ops.wipe(); T.define(ops)
base = json.load(open(T.BASELINE))["rows"]
jobs = T.jobs()
res = {}
for tag in T.PROTOS:
    rows = {}; raw = {}
    for kind, key, st, de in jobs:
        o = T.replay(ops, tag, st, de); rows[key] = T.row(o); raw[key] = o
    res[tag] = (rows, raw)
for tag in (11, 12):
    rows = res[tag][0]
    bad = [k for k in base if not T.rows_equal(rows[k], base[k])]
    print(f"proto {tag}: byte-identity vs baseline: {len(base)-len(bad)}/{len(base)} equal; first bad: {bad[:3]}")
codes = {1:"startF",2:"startAlpha",3:"startOther",4:"errorAtDTmin",5:"loadingNonPosH",6:"tension",7:"drift",8:"alphaOut",9:"maxSubsteps"}
for tag in T.PROTOS:
    rows, raw = res[tag]
    for kind in ("fan","ring","repro"):
        keys = [k for (kk,k,_,_) in jobs if kk == kind]
        ref = collections.Counter(codes.get(int(float.fromhex(rows[k][26])), "?") for k in keys if rows[k][0] != 0)
        sub = [raw[k]["sas"]["substeps"] for k in keys]
        rs = sum(raw[k]["sas"]["alphaInReseats"] for k in keys); rv = sum(raw[k]["sas"]["rejectedReversal"] for k in keys)
        extra = ""
        if raw[keys[0]]["sas"] and "hFloored" in raw[keys[0]]["sas"]:
            extra = f" floored {sum(raw[k]['sas']['hFloored'] for k in keys):.0f} softcapped {sum(raw[k]['sas']['hSoftCapped'] for k in keys):.0f} held {sum(raw[k]['sas']['reseatHeld'] for k in keys):.0f}"
        print(f"proto {tag} {T.PROTOS[tag]!s:60s} {kind:5s}: refused {sum(1 for k in keys if rows[k][0]!=0):3d}/{len(keys)} {dict(ref)} substeps med {statistics.median(sub):.0f} max {max(sub):.0f} reseats {rs:.0f} rejRev {rv:.0f}{extra}")
# oracle comparison on the fan
ora = json.load(open(T.ORACLE_FAN))["fan"]
for tag, ov in ((13, "T1B1S"), (16, "T1B1")):
    rels = []
    for kind, key, st, de in jobs:
        if kind != "fan": continue
        o = res[tag][1][key]; r = ora[key][ov]
        if o["rc"] != 0 or r["status"] != "ok": continue
        dc = [a - b for a, b in zip(o["sigma"], st["sigma"])]
        do = r["dsig"]
        nrm = lambda v: math.sqrt(v[0]**2+v[1]**2+v[2]**2+2*(v[3]**2+v[4]**2+v[5]**2))
        rels.append(nrm([a-b for a,b in zip(dc,do)]) / max(nrm(do), 1e-12))
    rels.sort()
    print(f"proto {tag} vs oracle {ov}: n={len(rels)} rel |dsig diff| median {statistics.median(rels):.2e} p95 {rels[int(0.95*len(rels))-1]:.2e} max {rels[-1]:.2e}")
