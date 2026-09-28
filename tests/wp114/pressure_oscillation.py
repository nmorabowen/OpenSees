"""WP-114 ask 3 evidence: intra-element pressure (I1) spread in the plastic
zone at the same punch settlement, pre- vs post-fix, psi=0 and psi=phi.
Locking under isochoric flow shows up as GP-to-GP I1 oscillation inside an
element (checkerboard); a relieved element keeps I1 smooth.

    python3.12 pressure_oscillation.py   (runs 4 subprocesses, s/B = 0.03)
"""
import json, os, subprocess, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
BINS = {'pre': os.path.join(HERE, 'prefix_bin'),
        'post': os.path.normpath(os.path.join(HERE, '..', '..', 'dist', 'bin'))}
OUT = os.path.join(HERE, 'out', 'osc')
os.makedirs(OUT, exist_ok=True)
for psi in ('0', 'phi'):
    for b in ('pre', 'post'):
        f = os.path.join(OUT, f'{b}_psi{psi}.json')
        if not os.path.exists(f):
            subprocess.run([sys.executable, os.path.join(HERE, 'punch_nonassoc.py'), '--bin', BINS[b],
                            '--elem', 'tri6_bbar', '--grid', '1', '--psi', psi, '--smax', '0.03',
                            '--out', f], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           env=dict(os.environ, LADRUNO_OPENSEES_QUIET='1'))
        r = json.load(open(f))
        by = {}
        for e, g, br, i1 in r['gp_all']:
            by.setdefault(e, []).append((br, i1))
        spread = [(max(v[1] for v in L) - min(v[1] for v in L)) / max(abs(np.mean([v[1] for v in L])), 1e-9)
                  for L in by.values() if any(v[0] == 1 for v in L)]
        sp = np.array(spread)
        print(f"{b:4s} psi={psi:3s} build={r['build'][:9]} s/B={r['s_over_B_last']:.4f} q={r['q_last']:.1f} "
              f"plastic elems={len(sp)} I1 spread/|mean| median={np.median(sp):.3f} p90={np.percentile(sp, 90):.3f} "
              f"cutoff+apex={r['n_cutoff'] + r['n_apex']}")
