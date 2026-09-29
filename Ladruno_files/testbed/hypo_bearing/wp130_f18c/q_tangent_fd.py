"""WP-130 / F18(c): is IntScheme 2's TanType-2 tangent the derivative of its own
return map?  (Why the bearing leg's global Newton diverges from iteration 1.)

3D cube (WP-127 byte-id deck, LadrunoSANISAND, IntScheme 2, TanType 2), 20
plastic steps; at committed step k the next step's strain increment de is
replayed from the committed state (ladrunoSANISANDReplay, exact to 1e-9 per
WP-127) with de +- h*e_j, giving the finite-difference algorithmic tangent
D_fd = d sigma_{k+1} / d eps_{k+1}. It is compared with the `tangent` response
the element reads after the analysis' own step k+1 (mCep_Consistent from
NewtonSol, TanType 2). Engineering shear strain; tension-positive stress.

Only rows whose replay reproduces the analysis' own k+1 stress (replay-vs-analysis
~1e-15) are a valid comparison; later rows of this driver drift (1e-1) and are
printed for completeness only.

Also the same for IntScheme 1 (ModifiedEuler's chained TanType-2 tangent) as a
control. Prints relative Frobenius errors ||D - D_fd|| / ||D_fd|| per scheme.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "Ladruno_scripts"))
sys.path.insert(0, os.path.join(ROOT, "tests"))
import sanisand_replay as sr  # noqa: E402
import wp127_sanisand_byteid as b127  # noqa: E402
from _testbed import ops  # noqa: E402


def grab():
    g = lambda name: list(ops.eleResponse(1, "material", 1, name))
    return dict(sig=g("stress"), eps=g("strain"), alpha=g("alpha"),
                ain=g("alpha_in"), z=g("fabric"), e=g("state")[24], tan=g("tangent"))


def check(scheme, extra=(), nsteps=20, h=1e-8, lat=0.5, e_ax=5e-3):
    incs = b127._iso_dev(40, e_ax, lat)
    b127._build_3d("LadrunoSANISAND", b127._PARAMS, (scheme, 2, 1, 1e-7, 1e-7) + tuple(extra),
                   10, 3.0e-6, incs)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(10):
        assert ops.analyze(1) == 0
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    out = []
    prev = grab()
    for step in range(nsteps):  # every 5th step is checked; ONE extra analyze per check
        assert ops.analyze(1) == 0
        k = grab()
        if step >= 2:
            k1_eps = None
        # replay the NEXT step from k: need k+1 strain -> take one more analyze
        if step % 5 == 4:
            assert ops.analyze(1) == 0
            k1 = grab()
            de = [a - b for a, b in zip(k1["eps"], k["eps"])]
            dnorm = math.sqrt(sum((a - b) ** 2 for a, b in zip(k["eps"][:3], prev["eps"][:3]))
                              + 0.5 * sum((a - b) ** 2 for a, b in zip(k["eps"][3:], prev["eps"][3:])))
            base = sr.replay(ops, 1, k["sig"], k["alpha"], k["ain"], k["z"], k["e"], de,
                             "tensionPositive", prev_incr_norm=dnorm)
            D = [[0.0] * 6 for _ in range(6)]
            for j in range(6):
                dp = list(de); dp[j] += h
                dm = list(de); dm[j] -= h
                rp = sr.replay(ops, 1, k["sig"], k["alpha"], k["ain"], k["z"], k["e"], dp,
                               "tensionPositive", prev_incr_norm=dnorm)
                rm = sr.replay(ops, 1, k["sig"], k["alpha"], k["ain"], k["z"], k["e"], dm,
                               "tensionPositive", prev_incr_norm=dnorm)
                for i in range(6):
                    D[i][j] = (rp["sigma"][i] - rm["sigma"][i]) / (2 * h)
            T = [[k1["tan"][6 * i + j] for j in range(6)] for i in range(6)]
            nfd = math.sqrt(sum(D[i][j] ** 2 for i in range(6) for j in range(6)))
            err = math.sqrt(sum((T[i][j] - D[i][j]) ** 2 for i in range(6) for j in range(6)))
            # the `tangent` response comes out with the OPPOSITE sign to a tension-positive
            # d sigma / d eps (measured: every entry); compare -T as well
            errT = math.sqrt(sum((-T[i][j] - D[i][j]) ** 2 for i in range(6) for j in range(6)))
            asym = math.sqrt(sum((D[i][j] - D[j][i]) ** 2 for i in range(6) for j in range(6))) / nfd
            reps = max(abs(a - b) for a, b in zip(base["sigma"], k1["sig"])) / max(abs(x) for x in k1["sig"])
            if step == 4 and os.environ.get("FD_VERBOSE"):
                print("T  rows 0-2", [["%.4g" % T[i][j] for j in range(6)] for i in range(3)])
                print("FD rows 0-2", [["%.4g" % D[i][j] for j in range(6)] for i in range(3)])
            out.append((step, err / nfd, errT / nfd, asym, reps, base["stats"]["cppmNewtonFail"], base["path"]))
            k = k1
        prev = k
    return out


for scheme, extra in [(2, ()), (1, ())]:
    for row in check(scheme, extra, nsteps=int(os.environ.get('FD_STEPS', 20))):
        print(f"IntScheme {scheme} step {row[0]:2d}: ||T-Dfd||/||Dfd|| = {row[1]:.3e} "
              f"||-T-Dfd||/||Dfd|| = {row[2]:.3e}, asym(Dfd) {row[3]:.2e}, replay-vs-analysis {row[4]:.1e}, "
              f"cppmNewtonFail {int(row[5])}, path {row[6]}", flush=True)
