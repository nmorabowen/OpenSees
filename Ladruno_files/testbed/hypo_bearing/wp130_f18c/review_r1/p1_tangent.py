import sys
exec(open(r"C:\Users\nmora\AppData\Local\Temp\claude\C--Users-nmora-Github-OpenSees-Compile-OpenSees--claude-worktrees-tims-implementation-review-3733c6\e034e494-142a-4e05-862e-9d267344b7ad\scratchpad\r1\h.py").read())

s3 = 1 / math.sqrt(3); s2 = 1 / math.sqrt(2)
DIRS = {"vol": [s3, s3, s3, 0, 0, 0], "dev12": [s2, -s2, 0, 0, 0, 0],
        "dev_ax": [-0.5 * s2 * 0 - 1 / math.sqrt(6), -1 / math.sqrt(6), 2 / math.sqrt(6), 0, 0, 0],
        "g12": [0, 0, 0, 1, 0, 0], "g23": [0, 0, 0, 0, 1, 0]}


def iso_dev(n, e_ax, lat):
    return [(-lat * e_ax / n, e_ax / n)] * n


def run_case(label, tolR, incs, e_conf, check_steps, tolF=1e-7, extra=(), h=1e-8):
    build_3d("LadrunoSANISAND", PARAMS, (2, 2, 1, tolF, tolR) + tuple(extra), 10, e_conf, incs)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(10):
        assert ops.analyze(1) == 0
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    hist = [grab()]
    for s in range(max(check_steps) + 1):
        rc = ops.analyze(1)
        assert rc == 0, (label, s, rc)
        hist.append(grab())
    for s in check_steps:
        prev, k, k1 = hist[s - 1], hist[s], hist[s + 1]
        names = list(DIRS)
        cols, Tcols, rep, base = fd_from(k, k1, prev, h=h, dirs=[DIRS[n] for n in names])
        p = -sum(k1["sig"][:3]) / 3
        errs = {n: relerr([Tc], [c]) for n, Tc, c in zip(names, Tcols, cols)}
        errs_m = {n: relerr([[-x for x in Tc]], [c]) for n, Tc, c in zip(names, Tcols, cols)}
        print(f"{label} tolR={tolR:g} step {s}: p={p:.3f} rep={rep:.1e} path={base['path']} "
              f"nfail={int(base['stats']['cppmNewtonFail'])} halv={int(base['stats']['cppmHalvings'])} expl={int(base['stats']['cppmExplicitFail'])} "
              + " ".join(f"{n}:{errs[n]:.2e}" for n in names), flush=True)


which = sys.argv[1] if len(sys.argv) > 1 else "a"
if which == "a":
    # the author's state (step 5 of iso_dev(40,5e-3,0.5), e_conf 3e-6)
    for tolR in (1e-7, 1e-10, 1e-12):
        run_case("author", tolR, iso_dev(40, 5e-3, 0.5), 3e-6, [5, 10, 20])
elif which == "b":
    # higher confinement, near-failure (large deviatoric strain)
    for tolR in (1e-7, 1e-11):
        run_case("conf1e-4", tolR, iso_dev(60, 3e-2, 0.5), 1e-4, [5, 30, 55])
elif which == "c":
    # reversal: load then reverse
    e = 2e-3
    incs = [(-0.5 * e / 10, e / 10)] * 10 + [(0.5 * e / 10, -e / 10)] * 15
    for tolR in (1e-7, 1e-11):
        run_case("reversal", tolR, incs, 1e-4, [9, 10, 11, 12, 20])
if which == "d":
    # high p (sigmoid inactive), Presidual 0 vs default
    for extra in ((), ("-Presidual", 0.0)):
        for tolR in (1e-7, 1e-12):
            run_case(f"hiP{extra}", tolR, iso_dev(40, 1e-2, 0.5), 3e-4, [3, 15, 35], extra=extra)
if which == "e":
    # low p (sigmoid active), Presidual 0
    for tolR in (1e-12,):
        run_case("lowP_pr0", tolR, iso_dev(40, 5e-3, 0.5), 3e-6, [5, 10], extra=("-Presidual", 0.0))
