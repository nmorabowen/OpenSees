"""WP-216 bench: per-mode RSA readback, today's loop against responseSpectrumPerMode.

    python tests/bench_wp216_rsa_per_mode.py --nst 10 --h 0.4 --modes 25
    python tests/bench_wp216_rsa_per_mode.py --nst 45 --h 0.5 --modes 60
    python tests/bench_wp216_rsa_per_mode.py --nst 45 --h 0.5 --modes 60 --skip-old --linear

Synthetic tower from tests/_testbed/rsa_tower.py (no project data). Today's
loop is the one the tower project ran: per direction and mode
`responseSpectrumAnalysis -mode k -increment`, `reactions`, `nodeReaction` on
the supports, `eleForce` on the pier shells, inertia sums in Python (from the
IncrDisp node response: `nodeDisp` is the committed state). Prints
seconds per (mode, direction) for every phase, then the new command.
--linear (WP-223) also times `-linear` against the default call and prints
the largest difference per output group, relative to its largest component.
Not collected by pytest (no test_ prefix).
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
from _testbed import ops                         # noqa: E402
from _testbed.rsa_tower import build_tower, spectrum   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nst', type=int, default=10)
    ap.add_argument('--h', type=float, default=0.4)
    ap.add_argument('--modes', type=int, default=25)
    ap.add_argument('--system', default='Pardiso')
    ap.add_argument('--loop-modes', type=int, default=0,
                    help='time the old loop on only this many modes (0 = all)')
    ap.add_argument('--skip-old', action='store_true')
    ap.add_argument('--breakdown', action='store_true',
                    help='also time the new command on reactions / piers / inertia alone')
    ap.add_argument('--linear', action='store_true',
                    help='WP-223: also time -linear against the default call')
    ap.add_argument('--preload', action='store_true',
                    help='gravity on every floor node, then loadConst (before eigen)')
    ap.add_argument('--transf', default='Linear')
    ap.add_argument('--pdelta-shells', action='store_true')
    ap.add_argument('--all-piers', action='store_true',
                    help='request the pier group of every storey, not storeys 1 and nst/2')
    a = ap.parse_args()

    try:
        print('build', ops.ladrunoBuild())
    except Exception:
        pass
    t0 = time.perf_counter()
    t = build_tower(a.nst, a.h, system=a.system, transf=a.transf,
                    shell_args=('-pDelta',) if a.pdelta_shells else ())
    if a.preload:
        ops.timeSeries('Constant', 1)
        ops.pattern('Plain', 1, 1)
        for s in t['storey_nodes']:
            for n in t['storey_nodes'][s]:
                if abs(ops.nodeCoord(n)[2] - 3.0 * s) < 1e-9:
                    ops.load(n, 0.0, 0.0, -2.0e4, 0.0, 0.0, 0.0)
        ops.integrator('LoadControl', 1.0)
        ops.analysis('Static')
        assert ops.analyze(1) == 0
        ops.loadConst('-time', 0.0)
    t1 = time.perf_counter()
    print(f"tower nst={a.nst} h={a.h}: {t['n_nodes']} nodes, {t['n_ele']} elements, "
          f"{len(t['supports'])} supports; build {t1 - t0:.1f} s")
    lam = ops.eigen('-genBandArpack', a.modes)
    t2 = time.perf_counter()
    ops.modalProperties('-unorm')
    t3 = time.perf_counter()
    print(f"eigen {a.modes}: {t2 - t1:.1f} s, modalProperties {t3 - t2:.2f} s, "
          f"T1 = {2 * np.pi / np.sqrt(lam[0]):.3f} s")

    tn, sa = spectrum()
    spec = ['-Tn', *tn, '-Sa', *sa]
    piers = t['piers']
    in_s = [1, max(1, a.nst // 2)]                  # inertia storeys
    pier_s = sorted(piers) if a.all_piers else in_s
    sup = t['supports']
    st = t['storey_nodes']
    dirs = (1, 2)
    nm = a.loop_modes or a.modes

    if not a.skip_old:
        tim = dict(rsa=0.0, reactions=0.0, nodeReaction=0.0, eleForce=0.0, inertia=0.0)
        mass = {n: np.array(ops.nodeMass(n)) for s in in_s for n in st[s]}
        for d in dirs:
            for k in range(1, nm + 1):
                c0 = time.perf_counter()
                ops.responseSpectrumAnalysis(d, *spec, '-mode', k, '-increment')
                c1 = time.perf_counter()
                ops.reactions()
                c2 = time.perf_counter()
                R = np.sum([ops.nodeReaction(n) for n in sup], axis=0)
                c3 = time.perf_counter()
                for s in pier_s:
                    for e in piers[s][0]:
                        ops.eleForce(e)
                c4 = time.perf_counter()
                for s in in_s:
                    f = np.zeros(6)
                    for n in st[s]:
                        # trial - committed (nodeDisp is the committed state)
                        du = [ops.nodeResponse(n, i, 4) for i in range(1, 7)]
                        f += lam[k - 1] * mass[n] * np.array(du)
                c5 = time.perf_counter()
                tim['rsa'] += c1 - c0
                tim['reactions'] += c2 - c1
                tim['nodeReaction'] += c3 - c2
                tim['eleForce'] += c4 - c3
                tim['inertia'] += c5 - c4
        ops.responseSpectrumAnalysis('-revert')
        n = len(dirs) * nm
        tot = sum(tim.values())
        print(f"OLD loop, {n} (mode, dir) pairs: {tot:.2f} s = {tot / n * 1e3:.1f} ms per pair")
        for key, v in tim.items():
            print(f"    {key:13s} {v / n * 1e3:8.1f} ms/pair  ({100 * v / tot:5.1f} %)")

    if hasattr(ops, 'responseSpectrumPerMode'):
        rea = ['-reactions']
        ele, ine = [], []
        for s in pier_s:
            e, c = piers[s]
            ele += ['-eleGroup', f'pier{s}', '-ele', *e, '-cut', *c]
        for s in in_s:
            ine += ['-inertia', f'storey{s}', *st[s], '-about', 15.0, 15.0, 0.0]
        n = len(dirs) * a.modes
        runs = [('all', rea + ele + ine)]
        if a.breakdown:
            runs += [('reactions only', rea), ('reactions + piers', rea + ele),
                     ('inertia only', ine)]
        if a.linear:
            runs += [('all -linear', rea + ele + ine + ['-linear'])]
        outs = {}
        for label, extra in runs:
            c0 = time.perf_counter()
            out = ops.responseSpectrumPerMode(*spec, '-dir', *dirs, *extra)
            c1 = time.perf_counter()
            outs[label] = out
            print(f"NEW responseSpectrumPerMode [{label}], {n} (mode, dir) pairs: {c1 - c0:.3f} s = "
                  f"{(c1 - c0) / n * 1e3:.2f} ms per pair; updated {int(out['numUpdated'][0])} "
                  f"of {t['n_ele']} elements")
        if a.linear:
            on, off = outs['all -linear'], outs['all']
            print(f"-linear check modes {[int(v) for v in on['linearCheckModes']]}, "
                  f"check error {on['linearCheckError'][0]:.2e}")
            nd, nm, nc = (int(v) for v in off['shape'])
            for key in off:
                if ':' in key or key == 'reactions':
                    x = np.array(off[key]).reshape(nd, nm, nc)
                    y = np.array(on[key]).reshape(nd, nm, nc)
                    # against the key's largest value of the kind over every
                    # (dir, mode): per (dir, mode) a near-zero Gamma compares
                    # round-off with round-off
                    e = 0.0
                    for lo, hi in ((0, 3), (3, 6)):
                        sc = np.abs(x[..., lo:hi]).max()
                        if sc > 0.0:
                            e = max(e, float(np.abs(y[..., lo:hi] - x[..., lo:hi]).max() / sc))
                    print(f"    {key:14s} max |on - off| / max |off| per kind: {e:.2e}")


if __name__ == '__main__':
    main()
