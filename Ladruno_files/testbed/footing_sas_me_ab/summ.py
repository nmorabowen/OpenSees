"""Summarise pulled esmeralda legs: last steps, cumulatives, E_A vs E_B curve."""
import csv, os, re, sys
H = os.path.dirname(os.path.abspath(__file__))
LEGS = ["E_A", "E_B", "E_C2", "E_B16"]
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4


def rows(l):
    p = os.path.join(H, "runs", l, "steps.csv")
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


for l in LEGS:
    R = rows(l)
    if not R:
        print(l, "no rows"); continue
    last = R[-1]
    tot = float(last["wall_total_s"])
    fails = sum(int(r["fails_before"]) for r in R)
    caps = sum(int(r["cap_step"]) for r in R)
    dtmin = sum(int(r["dtmin_step"]) for r in R)
    rungs = {k: sum(1 for r in R if r["rung"] == k) for k in "012"}
    qmax = max(float(r["q_kPa"]) for r in R)
    print(f"== {l}: steps {last['step']} s/B {float(last['s_over_B']):.5f} q {float(last['q_kPa']):.2f} "
          f"(qmax {qmax:.2f}) push wall {tot/3600:.2f} h; failed attempts {fails}; "
          f"refusals/caps {caps}; dtmin {dtmin}; rungs N/L/K {rungs['0']}/{rungs['1']}/{rungs['2']}")
    for r in R[-N:]:
        print(f"   {int(r['step']):4d} s/B {float(r['s_over_B']):.5f} q {float(r['q_kPa']):8.2f} ds {float(r['ds_m']):.2e} "
              f"rung {'NLK'[int(r['rung'])]} it {int(r['iters']):3d} fails {r['fails_before']} "
              f"wall {float(r['wall_step_s']):8.1f}s sub {int(r['sub_step_total']):.2e} cap {r['cap_step']} "
              f"dtmin {r['dtmin_step']} rho_max {r['max_rho_alpha']} n>1 {r['n_rho_gt_1']} p'min {r['min_p']}")
    lg = os.path.join(H, "runs", l, "logs", "log.log")
    txt = open(lg, errors="replace").read()
    for m in re.findall(r".*(MODE = .*|DEAD-END candidate.*|exit .* after.*|Traceback.*)", txt)[-4:]:
        print("   LOG:", m[:200])


def interp(R, s):
    xs = [float(r["s_over_B"]) for r in R]; qs = [float(r["q_kPa"]) for r in R]
    for i in range(1, len(xs)):
        if xs[i] >= s:
            t = (s - xs[i-1]) / (xs[i] - xs[i-1]); return qs[i-1] + t * (qs[i] - qs[i-1])
    return None


A, B = rows("E_A"), rows("E_B")
for ref, other, name in ((A, B, "E_B vs E_A"), (A, rows("E_C2"), "E_C2 vs E_A"), (B, rows("E_C2"), "E_C2 vs E_B"),
                         (rows("E_C2"), rows("E_B16"), "E_B16 vs E_C2")):
    if not ref or not other: continue
    smax = min(float(ref[-1]["s_over_B"]), float(other[-1]["s_over_B"]))
    worst = (0, 0, 0)
    for r in other:
        s = float(r["s_over_B"])
        if s > smax or s < float(ref[0]["s_over_B"]): continue
        qr = interp(ref, s)
        if qr is None: continue
        d = float(r["q_kPa"]) - qr
        if abs(d) > abs(worst[1]): worst = (s, d, d / qr)
    print(f"{name}: common s/B <= {smax:.5f}; worst dq {worst[1]:+.3f} kPa ({100*worst[2]:+.3f} %) at s/B {worst[0]:.5f}")
