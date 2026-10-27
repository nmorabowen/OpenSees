"""WP-216 bench: per-mode RSA readback, today's loop against responseSpectrumPerMode.

    python tests/bench_wp216_rsa_per_mode.py --nst 10 --h 0.4 --modes 25
    python tests/bench_wp216_rsa_per_mode.py --nst 45 --h 0.5 --modes 60

Synthetic tower from tests/_testbed/rsa_tower.py (no project data). Today's
loop is the one the tower project ran: per direction and mode
`responseSpectrumAnalysis -mode k -increment`, `reactions`, `nodeReaction` on
the supports, `eleForce` on the pier shells, inertia sums in Python. Prints
seconds per (mode, direction) for every phase, then the new command.
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
    a = ap.parse_args()

    try:
        print('build', ops.ladrunoBuild())
    except Exception:
        pass
    t0 = time.perf_counter()
    t = build_tower(a.nst, a.h, system=a.system)
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
    pier_s = [1, max(1, a.nst // 2)]
    sup = t['supports']
    st = t['storey_nodes']
    dirs = (1, 2)
    nm = a.loop_modes or a.modes

    if not a.skip_old:
        tim = dict(rsa=0.0, reactions=0.0, nodeReaction=0.0, eleForce=0.0, inertia=0.0)
        mass = {n: np.array(ops.nodeMass(n)) for s in pier_s for n in st[s]}
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
                for s in pier_s:
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
        args = [*spec, '-dir', *dirs, '-reactions']
        for s in pier_s:
            e, c = piers[s]
            args += ['-eleGroup', f'pier{s}', '-ele', *e, '-cut', *c]
            args += ['-inertia', f'storey{s}', *st[s], '-about', 15.0, 15.0, 0.0]
        c0 = time.perf_counter()
        out = ops.responseSpectrumPerMode(*args)
        c1 = time.perf_counter()
        n = len(dirs) * a.modes
        print(f"NEW responseSpectrumPerMode, {n} (mode, dir) pairs: {c1 - c0:.3f} s = "
              f"{(c1 - c0) / n * 1e3:.2f} ms per pair; updated {int(out['numUpdated'][0])} "
              f"of {t['n_ele']} elements")


if __name__ == '__main__':
    main()
