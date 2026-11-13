"""Siddiquee legs at the test's stages: python k5_sid_cmp.py <leg> ...  (run in ~/wp144/deck/runs)"""
import csv, json, sys
import numpy as np
print("test, 1 g: N 50.7 / 78.2 / 105.0 at 1.99 / 3.20 / 4.82 %, peak 110.6 at 5.22 %; centrifuge peak 159.7 at 12.0 %")
for d in sys.argv[1:]:
    try:
        R = list(csv.DictReader(open(d + "/steps.csv")))
    except Exception as e:
        print(d, "missing", e); continue
    try: sm = json.load(open(d + "/summary.json"))
    except Exception: sm = {}
    if not R: print(d, "no steps"); continue
    s = np.array([float(r["s_over_B"]) for r in R]) * 100; N = np.array([float(r["q_kPa"]) for r in R]) / (15.6 * 0.5 / 2)
    st = [float(np.interp(x, s, N)) if s[-1] >= x else float("nan") for x in (1.99, 3.20, 4.82)]
    print("%-22s N %5.1f %5.1f %5.1f (x %.2f %.2f %.2f) | max %6.1f at %4.1f %% (x %.2f / %.2f) | end %.1f %% mode %s fails %s refused %s wall %d" % (
        d, *st, st[0] / 50.7, st[1] / 78.2, st[2] / 105.0, N.max(), s[N.argmax()], N.max() / 110.6, N.max() / 159.7, s[-1], sm.get("mode"),
        sm.get("failed_attempts"), (sm.get("norsand") or {}).get("refused_trials"), sm.get("total_wall_s") or 0))
