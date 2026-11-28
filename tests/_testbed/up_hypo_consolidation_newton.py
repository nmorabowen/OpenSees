"""WP-266 (D-UP1) -- large-strain consolidation under LadrunoUP -geom hypo:
Newton iterations per step with the shipped (dominant-term) tangent, and a
byte-exact trajectory digest.

WHAT IT RUNS
    An oedometric column of H8 LadrunoUP elements (-geom hypo, -kozenyCarman,
    ElasticIsotropic skeleton), drained at the top, loaded by a sudden surcharge
    q = --qfrac * E: the undrained response first, then consolidation (defaults
    q/E = 0.2, 40 steps of dt = 0.2: about 10 % mean vertical strain, more at the
    drained top). Each variant is the same deck under a different integrator or
    element-argument set:

      newmark_battery   Newmark(0.6, 0.3025), the -geom hypo battery's integrator.
                        The integrator scales the element's p-row rate block
                        [Q_cur', S + aH~] by c2 = gamma/(beta dt) = 1.98/dt; the
                        incremental residual differentiates to 1/dt. The WP-266
                        oracle measures this as an O(1) error on the p-rows.
      newmark_matched   Newmark(0.5, 0.5): c2 = 1/dt, the rate block exact; what
                        is left is the dropped Jacobian blocks only.
      <name>=<args>     --variant adds one more: an integrator "gamma,beta" and
                        extra element arguments, e.g. the wave-2 P4 full-Jacobian
                        flag once it exists:
                        --variant "full=0.5,0.5,-tangent,full"

    Per variant it prints, on stdout: Newton iterations per step (ops.testIter),
    total iterations, failed steps, wall time, the final settlement and base
    excess pressure. The ladrunoBuild() stamp and the engine path go first.

DIGEST MODE (--digest)
    Prints ONLY the physics of the newmark_battery run, one value per line as
    float.hex() (top settlement and base p after every step, then every nodal
    displacement and pressure at the end). The default tangent path is unchanged
    by WP-266, so the digest of a pre-WP-266 build and of the WP-266 build must be
    IDENTICAL AS BYTES. Build stamps go to stderr.

USAGE (the orchestrator, after building; from tests/)
    py -3.12 -S _testbed/up_hypo_consolidation_newton.py
    py -3.12 -S _testbed/up_hypo_consolidation_newton.py --digest > wp266.digest.txt
    py -3.12 -S _testbed/up_hypo_consolidation_newton.py --digest --dist <pre-WP dist/bin> > pre.digest.txt
    fc /b wp266.digest.txt pre.digest.txt

    exit 0: every variant ran (a failed step is a RESULT, recorded, not an error)
    exit 2: no engine to bind

DRY RUN (the WP-266 worker, 2026-11-28, on the PRE-WP-266 engine
C:/Users/nmora/Documents/Github/OpenSees_ladruno/dist/bin, ladrunoBuild
f0b9a1df699f9634fd89c5942d69d90fe5ad27c4 -- the default tangent path WP-266
leaves unchanged; the orchestrator re-runs on the WP-266 build):
    defaults           newmark_battery 1045 iterations / 40 steps (26.1 per step),
                       newmark_matched  385 (9.6 per step), no failures; settlement
                       -0.3839 vs -0.3839.
    --qfrac 0.25 --dt 0.05 (the hard case, --print 1): newmark_battery contracts
                       the residual by 0.49 per iteration (= 1 - beta/gamma, the c2
                       factor), newmark_matched by 0.22 (the dropped blocks); at step
                       4 / step 6 Newton DIVERGES (the residual grows from 5e-4 to
                       3e3) -- the stall the full Jacobian is meant to remove.
"""
import argparse
import os
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve()
_TESTS = _HERE.parents[1]
_ROOT = _HERE.parents[2]

os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")


def _bind(dist):
    dist = Path(dist)
    if not ((dist / "opensees.pyd").is_file() or (dist / "opensees.so").is_file()):
        sys.stderr.write("up_hypo_consolidation_newton: FATAL - no engine at %s\n" % dist)
        raise SystemExit(2)
    if str(_TESTS) not in sys.path:
        sys.path.insert(0, str(_TESTS))
    from _engine import bind_worktree_engine
    return bind_worktree_engine(str(dist))


# --------------------------------------------------------------------------
# the deck
# --------------------------------------------------------------------------
E = 1.0e4
NU = 0.3
RHO_SAT = 2.0
KF = 2.2e6
PORO = 0.4
PERM = 1.0e-4


def _column(ops, nz, lz):
    tags, t = {}, 1
    for k in range(nz + 1):
        z = lz * k / nz
        for (x, y) in ((0, 0), (1, 0), (1, 1), (0, 1)):
            ops.node(t, float(x), float(y), z)
            tags[(k, (x, y))] = t
            t += 1
    conns = []
    for k in range(nz):
        conns.append([tags[(k, (0, 0))], tags[(k, (1, 0))], tags[(k, (1, 1))],
                      tags[(k, (0, 1))], tags[(k + 1, (0, 0))], tags[(k + 1, (1, 0))],
                      tags[(k + 1, (1, 1))], tags[(k + 1, (0, 1))]])
    return tags, conns


def run(ops, gamma, beta, extra, nz, lz, nsteps, dt, qfrac, tol, maxit, test="NormUnbalance", pflag=0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    tags, conns = _column(ops, nz, lz)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO_SAT)
    for e, conn in enumerate(conns, start=1):
        ops.element("LadrunoUP", e, *conn, 1, "-Kf", KF, "-poro", PORO, "-rhoF", 1.0,
                    "-perm", PERM, PERM, PERM, "-geom", "hypo", "-kozenyCarman", *extra)
    top = nz
    for (k, xy), t in tags.items():
        ops.fix(t, 1, 1, 1 if k == 0 else 0, 1 if k == top else 0)
    q = qfrac * E
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for xy in ((0, 0), (1, 0), (1, 1), (0, 1)):
        ops.load(tags[(top, xy)], 0.0, 0.0, -q / 4.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test(test, tol, maxit, pflag)
    ops.algorithm("Newton")
    ops.integrator("Newmark", gamma, beta)
    ops.analysis("Transient")

    its, traj, failed = [], [], 0
    t0 = time.perf_counter()
    for _ in range(nsteps):
        rc = ops.analyze(1, dt)
        its.append(ops.testIter() if rc == 0 else -ops.testIter())
        if rc != 0:
            failed += 1
            break
        traj.append(ops.nodeDisp(tags[(top, (0, 0))], 3))
        traj.append(ops.nodeDisp(tags[(0, (0, 0))], 4))
    wall = time.perf_counter() - t0
    final = []
    for t in sorted(tags.values()):
        final.extend(ops.nodeDisp(t))
    return its, traj, final, failed, wall


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dist", default=str(_ROOT / "dist" / "bin"),
                    help="engine directory (default: THIS worktree's dist/bin)")
    ap.add_argument("--nz", type=int, default=8)
    ap.add_argument("--lz", type=float, default=4.0)
    ap.add_argument("--nsteps", type=int, default=40)
    ap.add_argument("--dt", type=float, default=0.2)
    ap.add_argument("--qfrac", type=float, default=0.2,
                    help="surcharge q / E (the hard case: 0.25 with --dt 0.05)")
    ap.add_argument("--tol", type=float, default=1.0e-8)
    ap.add_argument("--maxit", type=int, default=100)
    ap.add_argument("--test", default="NormUnbalance", help="OpenSees convergence test")
    ap.add_argument("--print", dest="pflag", type=int, default=0, help="test printFlag (2: every iterate)")
    ap.add_argument("--variant", action="append", default=[],
                    help='extra variant "name=gamma,beta[,elementArg,...]"')
    ap.add_argument("--digest", action="store_true",
                    help="print only the newmark_battery physics as float.hex()")
    a = ap.parse_args(argv)

    ops = _bind(a.dist)
    build = str(ops.ladrunoBuild()) if hasattr(ops, "ladrunoBuild") else "unknown"
    deck = (a.nz, a.lz, a.nsteps, a.dt, a.qfrac, a.tol, a.maxit, a.test, a.pflag)

    if a.digest:
        sys.stderr.write("ladrunoBuild %s  engine %s\n" % (build, a.dist))
        its, traj, final, failed, _ = run(ops, 0.6, 0.3025, [], *deck)
        for v in its:
            print("iter %d" % v)
        for v in traj:
            print(float(v).hex())
        for v in final:
            print(float(v).hex())
        print("failed %d" % failed)
        return 0

    variants = [("newmark_battery", 0.6, 0.3025, []), ("newmark_matched", 0.5, 0.5, [])]
    for spec in a.variant:
        name, rest = spec.split("=", 1)
        parts = rest.split(",")
        extra = []
        for x in parts[2:]:
            try:
                extra.append(float(x))
            except ValueError:
                extra.append(x)
        variants.append((name, float(parts[0]), float(parts[1]), extra))

    print("ladrunoBuild %s" % build)
    print("engine %s" % a.dist)
    print("deck nz=%d lz=%g nsteps=%d dt=%g q/E=%g tol=%g maxit=%d test=%s (H8 hypo, -kozenyCarman, "
          "E=%g nu=%g k=%g n0=%g Kf=%g)" % (deck[:8] + (E, NU, PERM, PORO, KF)))
    for name, g, b, extra in variants:
        its, traj, final, failed, wall = run(ops, g, b, extra, *deck)
        total = sum(abs(i) for i in its)
        settle = traj[-2] if len(traj) >= 2 else float("nan")
        pbase = traj[-1] if traj else float("nan")
        print("VARIANT %s gamma=%g beta=%g c2*dt=%.4f extra=%s" % (name, g, b, g / b, extra))
        print("  iters/step %s" % " ".join(str(i) for i in its))
        print("  total %d  steps %d  failed %d  mean %.2f  wall %.3f s  settlement %.6e  base p %.6e"
              % (total, len(its), failed, total / max(len(its), 1), wall, settle, pbase))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
