"""K5 elastic-stiffness test: compare the scaled legs with the reference legs. Run in ~/wp144/deck/runs."""
import csv
import json
import re
for d in ["K5r4_kimura_V85.6_b8_F", "K5el_kimura_V85.6_b8_G3", "K5r4_kimura_V75.1_b8_F", "K5el_kimura_V75.1_b8_G0.33",
          "K5el_kimura_V75.1_b8_G3", "K5el_kimura_V75.1_b8_G10", "K5r4_kimura_V64.6_b8_F", "K5el_kimura_V64.6_b8_G3"]:
    try:
        R = list(csv.DictReader(open(d + "/steps.csv")))
        sm = json.load(open(d + "/summary.json"))
    except Exception as e:
        print(d, "missing", e)
        continue
    S = [float(r["s_over_B"]) for r in R]
    Q = [float(r["q_kPa"]) for r in R]
    at = lambda x: Q[min(range(len(S)), key=lambda i: abs(S[i] - x))]
    qm = max(Q)
    im = Q.index(qm)
    half = next(S[i] for i in range(len(Q)) if Q[i] >= 0.5 * qm)
    m = re.search(r"'-k', ([0-9.]+)\], \['-g', ([0-9.]+)", open(d + "/logs/log.log").read())
    print("%-27s k %7.1f g %7.1f | peak %7.1f at %.4f | s/B at half peak %.4f | q at .005 .01 .02 .05 .10: %.0f %.0f %.0f %.0f %.0f | "
          "mode %s fails %s refused %s wall %d" % (d, float(m.group(1)), float(m.group(2)), qm, S[im], half, at(.005), at(.01),
                                                    at(.02), at(.05), at(.10), sm.get("mode"), sm.get("failed_attempts"),
                                                    (sm.get("norsand") or {}).get("refused_trials"), sm.get("total_wall_s") or 0))
