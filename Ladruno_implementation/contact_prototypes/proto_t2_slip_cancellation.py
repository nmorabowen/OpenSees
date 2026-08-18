import numpy as np
# Hypothesis: storing slip as a GLOBAL-FRAME VECTOR makes gT-gp suffer catastrophic
# cancellation, destroying the relative accuracy of the slip increment (s-sp).
th = np.linspace(-np.pi, np.pi, 200001)[1:]
c, sn = np.cos(th), np.sin(th)
tx, ty = -sn, c
rng = np.random.default_rng(0)
sp = rng.uniform(-1, 1, tx.size) * 10.0 ** rng.uniform(-2, 2, tx.size)
ds = rng.uniform(-1, 1, tx.size) * 10.0 ** rng.uniform(-16, -8, tx.size)
s = sp + ds
# PATH A: slip lives as a global-frame vector, differenced componentwise
dA = (s * tx - sp * tx) * tx + (s * ty - sp * ty) * ty
# PATH B: slip lives as a scalar in the tangential frame
dB = s - sp
tru = dB  # fl(s-sp) is exact by Sterbenz for close operands; reference
relA = np.abs((dA - tru) / np.where(tru == 0, 1, tru))
relB = np.abs((dB - tru) / np.where(tru == 0, 1, tru))
print("slip increment, RELATIVE error vs the true difference:")
print("  PATH A (global-frame vector store): median %.3e   max %.3e   >100%% err in %.1f%% of cases"
      % (np.median(relA), relA.max(), 100 * np.mean(relA > 1.0)))
print("  PATH B (tangential scalar store)  : median %.3e   max %.3e"
      % (np.median(relB), relB.max()))
print()
print("  => cancellation scales with |ds|/|s|:")
ratio = np.abs(ds) / np.maximum(np.abs(sp), 1e-300)
for lo, hi in [(-16, -14), (-14, -12), (-12, -10), (-10, -8)]:
    m = (ratio > 10.0 ** lo) & (ratio <= 10.0 ** hi)
    if m.sum() > 50:
        print("     |ds|/|s| in [1e%d,1e%d]: PATH A median rel err %.2e  (n=%d)"
              % (lo, hi, np.median(relA[m]), m.sum()))
