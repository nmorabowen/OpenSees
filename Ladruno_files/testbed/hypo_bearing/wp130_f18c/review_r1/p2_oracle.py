import json, math, sys
R = r"C:\Users\nmora\AppData\Local\Temp\claude\C--Users-nmora-Github-OpenSees-Compile-OpenSees--claude-worktrees-tims-implementation-review-3733c6\e034e494-142a-4e05-862e-9d267344b7ad\scratchpad\r1"
sys.path.insert(0, R + r"\oracle")
from sanisand_reference import CAMPAIGN, Options, State, integrate
from sanisand_reference.model import bounding_report, v2t

rows = json.load(open(R + r"\p2_out_%s.json" % sys.argv[1]))
O = Options.uw(p_min=0.0101, p_residual=0.0, start_outside="plastic")
def dist(x, y):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(x, y)))
print("row  tag1err  tag2err  tag3err  tag4err | rho_end(ora,t1,t2,t3) | p(ora,t1,t2) | stats1 | status")
for i, row in enumerate(rows):
    r = row["res"]
    if not (r["2"]["stats"]["cppmGuessOk"] > 0 or r["3"]["stats"]["cppmGuessOk"] > 0):
        continue
    st = State.from_voigt(row["sig"], row["alpha"], row["z"], row["e"], row["ain"])
    try:
        res = integrate(st, row["de"], CAMPAIGN, O)
    except Exception as ex:
        print(i, "oracle error", ex); continue
    s = res.summary()
    ref = s["sigma"]
    sc = max(1.0, max(abs(x) for x in ref))
    errs = []
    rhos = [s["rho_alpha_end"]]
    for t in ("1", "2", "3", "4"):
        if r[t]["rc"] != 0:
            errs.append(float("nan"));
            if t in ("1", "2", "3"): rhos.append(float("nan"))
            continue
        errs.append(dist(r[t]["sigma"], ref) / sc)
        if t in ("1", "2", "3"):
            br = bounding_report(v2t(r[t]["sigma"]), v2t(r[t]["alpha"]), v2t(r[t]["z"]),
                                 row["e"] - (1 + 0.6944) * sum(row["de"][:3]), v2t(row["ain"]), CAMPAIGN, O)
            rhos.append(br["rho_alpha"])
    st1 = r["1"]["stats"]
    print(f"{i:3d} " + " ".join(f"{e:8.3f}" for e in errs) + " | " + " ".join(f"{x:6.3f}" for x in rhos)
          + f" | {s['p_end']:.2f} {r['1']['p']:.2f} {r['2']['p']:.2f} | nf{int(st1['cppmNewtonFail'])} h{int(st1['cppmHalvings'])} x{int(st1['cppmExplicitFail'])} | {s['status']}",
          flush=True)
