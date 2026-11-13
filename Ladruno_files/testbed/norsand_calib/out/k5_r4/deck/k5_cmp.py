"""Compare strip legs: python k5_cmp.py <leg> <leg> ...  (run in ~/wp144/deck/runs)"""
import csv, json, sys
for d in sys.argv[1:]:
    try:
        R = list(csv.DictReader(open(d + "/steps.csv"))); sm = json.load(open(d + "/summary.json"))
    except Exception as e:
        print(d, "missing", e); continue
    S = [float(r["s_over_B"]) for r in R]; Q = [float(r["q_kPa"]) for r in R]
    at = lambda x: Q[min(range(len(S)), key=lambda i: abs(S[i] - x))] if S[-1] >= x - 1e-3 else float("nan")
    qm = max(Q); im = Q.index(qm)
    print("%-27s peak %7.1f at %.4f | q at .005 .01 .02 .05 .10: %.0f %.0f %.0f %.0f %.0f | end %.3f mode %s fails %s refused %s wall %d" % (
        d, qm, S[im], at(.005), at(.01), at(.02), at(.05), at(.10), S[-1], sm.get("mode"), sm.get("failed_attempts"),
        (sm.get("norsand") or {}).get("refused_trials"), sm.get("total_wall_s") or 0))
