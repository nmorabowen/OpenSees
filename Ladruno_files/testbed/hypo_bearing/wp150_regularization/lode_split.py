"""WP-150: are the GPs that lose ellipticity on the EXTENSION side of the Lode interpolation?
The campaign c = 0.71 < 7/9 makes the DM04 g(theta, c) meridian non-convex near extension (WP-151 section 6.3).
cos3theta of n: +1 = triaxial compression, -1 = triaxial extension, 0 = plane-strain-like (the Lode angle is 0)."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import acoustic_vec as av

for f in sys.argv[1].split(","):
    o = av.run(f)
    bad = o["r"] <= 0
    c3 = o["c3"]
    tag = os.path.basename(os.path.dirname(os.path.dirname(f))) + "/" + os.path.basename(f)
    print(f"\n== {tag} s/B {o['sb']:.4f}: det<=0 {int(bad.sum())} of {len(c3)}")
    for name, m in (("det<=0", bad), ("elliptic", ~bad)):
        q = np.quantile(c3[m], [0.05, 0.25, 0.5, 0.75, 0.95])
        print(f"   {name:9s} cos3theta(n) quantiles 5/25/50/75/95 %: " + " ".join(f"{v:+.2f}" for v in q)
              + f"   share with cos3theta < -0.5 (extension side): {np.mean(c3[m] < -0.5):.2%}")
    edges = [-1.0, -0.5, 0.0, 0.5, 1.0001]
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (c3 >= lo) & (c3 < hi)
        if m.sum():
            print(f"   cos3theta in [{lo:+.1f},{hi:+.1f}): {int(m.sum()):5d} GPs, det<=0 fraction {np.mean(bad[m]):.2%}")
