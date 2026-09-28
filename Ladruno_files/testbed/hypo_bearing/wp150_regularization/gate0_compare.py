"""GATE 0 comparison table: C++ SAS-ME (R1 off / on) against the exact oracle (uw_model, the C++'s own equations), and
R1 on against R1 off. Metric: max |dq| over the path at matched axial strain, divided by q_max of the reference."""
import json, sys
import numpy as np

orc = json.load(open(sys.argv[1]))
cxx = json.load(open(sys.argv[2]))
ref_set = sys.argv[3] if len(sys.argv) > 3 else "uw_model"


def dq(a, b):
    ea, qa = np.array(a["eps_a"]), np.array(a["q"])
    eb, qb = np.array(b["eps_a"]), np.array(b["q"])
    grid = np.linspace(0, min(ea[-1], eb[-1]), 400)
    d = np.abs(np.interp(grid, ea, qa) - np.interp(grid, eb, qb))
    return float(d.max() / max(qa.max(), 1e-9))


print(f"build {cxx.get('build')}\n")
print(f"| test | oracle {ref_set} q_peak / q_end / e_end | C++ R1 off vs oracle max\|Δq\|/q_max | C++ R1 on vs off max\|Δq\|/q_max | R1 on vs off Δe_end |")
print("|---|---|---|---|---|")
worst_off, worst_on = 0.0, 0.0
for t, o in orc[ref_set].items():
    if t not in cxx["cxx_off"]:
        continue
    a, b = cxx["cxx_off"][t], cxx["cxx_on"][t]
    d1, d2 = dq(o, a), dq(a, b)
    worst_off, worst_on = max(worst_off, d1), max(worst_on, d2)
    print(f"| {t} | {max(o['q']):.1f} / {o['q'][-1]:.1f} / {o['e'][-1]:.4f} | {d1:.2e} | {d2:.2e} | {b['e'][-1]-a['e'][-1]:+.1e} |")
print(f"\nworst: C++ off vs oracle {worst_off:.2e}; R1 on vs off {worst_on:.2e}")
