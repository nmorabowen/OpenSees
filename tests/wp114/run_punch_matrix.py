"""WP-114: run the punch case matrix (one subprocess per case, so the pre- and
post-fix binaries never share a process) and print the report table.

    python3.12 run_punch_matrix.py            # runs the missing cases, prints table
    python3.12 run_punch_matrix.py --table    # table only, from out/punch_*.json
"""
import concurrent.futures as cf
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BINS = {'pre': os.path.join(HERE, 'prefix_bin'),
        'post': os.path.normpath(os.path.join(HERE, '..', '..', 'dist', 'bin'))}
OUT = os.path.join(HERE, 'out')

CASES = []
for psi in ('0', 'phi'):
    for b in ('pre', 'post'):
        CASES.append((b, 'tri6_bbar', 1, psi, 'bernstein'))
    CASES.append(('post', 'tri6_std', 1, psi, 'bernstein'))
    CASES.append(('post', 'quad', 1, psi, 'bernstein'))
    CASES.append(('post', 'quad', 2, psi, 'bernstein'))
    CASES.append(('post', 'tri6_bbar', 2, psi, 'bernstein'))
    # loading-artefact probe: Lagrange-weighted surcharge on Bernstein CPs
    CASES.append(('post', 'tri6_bbar', 1, psi, 'lagrange'))
    CASES.append(('pre', 'tri6_bbar', 1, psi, 'lagrange'))


def fname(c):
    return os.path.join(OUT, 'punch_%s_%s_g%d_psi%s_%s.json' % c)


def run(c):
    b, e, g, psi, load = c
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET='1')
    cmd = [sys.executable, os.path.join(HERE, 'punch_nonassoc.py'), '--bin', BINS[b],
           '--elem', e, '--grid', str(g), '--psi', psi, '--load', load, '--out', fname(c)]
    with open(fname(c).replace('.json', '.log'), 'w') as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env, timeout=3000)
    return c


def row(c):
    try:
        r = json.load(open(fname(c)))
    except FileNotFoundError:
        return '%-4s %-9s g%d psi=%-3s %-9s  (missing)' % c
    tag = '%-4s %-9s g%d psi=%-3s %-9s' % c
    if 'error' in r:
        return tag + '  STAGE-1 FAIL at t=%.2f; trial cutoff/apex GPs: %s' % (
            r['t_last_converged'], r['cutoff_apex_trial_xy_branch'])
    cen = r['census']
    rc = sorted(x['r_corner'] for x in cen)
    loc = ('r_corner min/med/max %.2f/%.2f/%.2f' % (rc[0], rc[len(rc) // 2], rc[-1])) if rc else ''
    return tag + ('  %s s/B=%.4f q_last=%7.1f q_peak=%7.1f@%.4f steps=%d cuts=%d fb=%d '
                  'cone=%d cutoff=%d apex=%d %s [%s %.0fs]') % (
        'WALL' if r['walled'] else 'ok  ', r['s_over_B_last'], r['q_last'], r['q_peak'],
        r['s_over_B_peak'], r['n_steps'], r['n_cuts'], r['n_alg_fallbacks'], r['n_cone'],
        r['n_cutoff'], r['n_apex'], loc, r['build'][:9], r['wall_s'])


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    if '--table' not in sys.argv:
        only = [a.split('=')[1] for a in sys.argv if a.startswith('--bins=')]
        todo = [c for c in CASES if not os.path.exists(fname(c))
                and (not only or c[0] in only[0].split(','))]
        with cf.ThreadPoolExecutor(max_workers=4) as ex:
            for c in ex.map(run, todo):
                print('done', c, flush=True)
    for c in CASES:
        print(row(c))
