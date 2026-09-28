"""The R1 copy with every toggle OFF must equal the committed WP-134 oracle."""
import sys, os, time, importlib
import numpy as np
W = sys.argv[1]; S = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(W, "Ladruno_scripts")); sys.path.insert(0, S)
import sanisand_reference as ref
import sanisand_r1 as r1
from sanisand_r1 import ring as r1ring
from sanisand_reference import ring as rring
CSV8 = os.path.join(W, "Ladruno_implementation/_tims_2d_model_requests_2026-09-25/ring_points_b8.csv")
rows = rring.load_ring_csv(CSV8)
cases = []
cases.append(("reproducer", rring.reproducer_state(), rring.REPRODUCER_DEPS))
for row in rows[:6]:
    for pn, de in rring.probes(1e-4).items():
        cases.append((f"b8 {row['element']}/{row['gp']} {pn}", rring.row_state(row), de))
opt_ref = rring.ring_variants()["uw_model"]
opt_r1 = r1ring.ring_variants()["uw_model"]
worst = 0.0
for name, st, de in cases:
    t0 = time.time()
    a = ref.integrate(st, de, ref.CAMPAIGN, opt_ref)
    st1 = r1.State(st.sigma.copy(), st.alpha.copy(), st.z.copy(), st.e, st.alpha_in.copy())
    b = r1.integrate(st1, de, r1.CAMPAIGN, opt_r1)
    d = max(np.max(np.abs(a.state.sigma - b.state.sigma)), np.max(np.abs(a.state.alpha - b.state.alpha)))
    worst = max(worst, d)
    print(f"{name:28s} {a.status:>14s} {b.status:>14s}  max|diff| {d:.1e}  ({time.time()-t0:.1f}s)")
print("WORST", worst)
