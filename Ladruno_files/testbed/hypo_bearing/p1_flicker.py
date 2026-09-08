"""ADR-95 P1: per-GP branch flicker between adjacent stations (0<->1 switches)."""
import sys, numpy as np
f = sys.argv[1]; z = np.load(f, allow_pickle=True)
sb = z['station_s_over_B']; n = len(sb)
B = np.stack([z[f'st{i:03d}_branch'] for i in range(n)])
F = np.stack([z[f'st{i:03d}_f1_trial'] for i in range(n)])
D = np.stack([z[f'st{i:03d}_detAmin'] for i in range(n)])
yielded = (B > 0).sum(1); sw = np.array([0] + [(B[i] != B[i-1]).sum() for i in range(1, n)])
ever = (B > 0).any(0).sum(); q = n // 4
core = (B[-q:] > 0).all(0).sum()
print(f'{f}: stations {n}, GPs {B.shape[1]}, ever-yielded {ever}, always-yielded-last-quarter {core}, switch rate mean {sw[1:].mean():.1f}/station')
print('  i    s/B       yielded switch  |f1|<1e-3  detA<0  f1max')
step = max(1, n // 18)
for i in sorted(set(list(range(0, n, step)) + [n-1])):
    print(f'{i:4d} {sb[i]:.6f} {yielded[i]:7d} {sw[i]:6d} {(np.abs(F[i])<1e-3).sum():9d} {(D[i]<=0).sum():7d} {F[i].max():8.4f}')
