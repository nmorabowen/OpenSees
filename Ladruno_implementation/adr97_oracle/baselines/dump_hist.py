"""ADR-94 wp/94c byte-identity gate -- FRESH SUBPROCESS per deck.

usage (cwd must be <worktree>/tests):
    PYTHONPATH=../dist/bin python3.12 dump_hist.py <out.json>      # driver
    PYTHONPATH=../dist/bin python3.12 dump_hist.py --one <deck> <tmp.json>
"""
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.abspath(os.getcwd()))


def _decks():
    """name -> (module_attr_spec) built lazily inside the child."""
    return [
        "tet/mc/BE/default", "tet/mc/BE/strict1", "tet/mc/BE/niter200",
        "tet/vm/BE/default", "tet/vm/BE/strict1",
        "tet/mctc/BE/default", "tet/mc_mctcfile/BE/default",
        "tet/hb/BE/default", "tet/dp/BE/default",
    ] + [
        "cube/vm/BE/%s/%s" % (t, leg)
        for t in ("Continuum", "Secant", "Elastic", "Numerical_Algorithmic_FirstOrder")
        for leg in ("plastic", "elastic")
    ] + [
        "cube/vm/BE/Continuum/soft", "cube/vm/BE/Continuum/hard0",
    ] + [
        "cube/vm/%s/Continuum" % m
        for m in ("Forward_Euler", "Forward_Euler_Subincrement",
                  "Modified_Euler_Error_Control", "Runge_Kutta_45_Error_Control")
    ]


def _build_one(name):
    """Return (build_callable, nsteps) for a deck name.  Child-process only."""
    from _testbed import ops  # noqa: F401
    import test_adr84_p2a_strict_convergence as P
    import test_asdplastic_mctc as M
    import test_adr94_hlist_numerics as N
    import test_adr94_hlist_hb as HB

    def tet(mat_fn):
        return (lambda: P._tet_build(mat_fn)), P.TET_NSTEPS

    def cube(mat_fn, load, nsteps):
        return (lambda: N._cube_build(mat_fn, load, nsteps)), nsteps

    table = {
        "tet/mc/BE/default": lambda: tet(lambda t: P.mat_mc(t)),
        "tet/mc/BE/strict1": lambda: tet(lambda t: P.mat_mc(t, strict=1)),
        "tet/mc/BE/niter200": lambda: tet(lambda t: P.mat_mc(t, niter=200)),
        "tet/vm/BE/default": lambda: tet(lambda t: P.mat_vm(t)),
        "tet/vm/BE/strict1": lambda: tet(lambda t: P.mat_vm(t, strict=1)),
        "tet/mctc/BE/default": lambda: tet(lambda t: M.mat_mctc(t)),
        "tet/mc_mctcfile/BE/default": lambda: tet(lambda t: M.mat_mc(t)),
        "tet/hb/BE/default": lambda: tet(lambda t: HB.mat_hb(t)),
        "tet/dp/BE/default": lambda: tet(lambda t: HB.mat_dp(t)),
        "cube/vm/BE/Continuum/soft":
            lambda: cube(lambda t: N.mat_vm(t, "Continuum", hiso=N.H_SOFT), N.P_PLASTIC, 6),
        "cube/vm/BE/Continuum/hard0":
            lambda: cube(lambda t: N.mat_vm(t, "Continuum", hiso=0.0), N.P_PLASTIC, 10),
    }
    for t in ("Continuum", "Secant", "Elastic", "Numerical_Algorithmic_FirstOrder"):
        table["cube/vm/BE/%s/plastic" % t] = (
            lambda g=t: cube(lambda t2, g=g: N.mat_vm(t2, g), N.P_PLASTIC, 10))
        table["cube/vm/BE/%s/elastic" % t] = (
            lambda g=t: cube(lambda t2, g=g: N.mat_vm(t2, g), -10.0, 10))
    for meth in ("Forward_Euler", "Forward_Euler_Subincrement",
                 "Modified_Euler_Error_Control", "Runge_Kutta_45_Error_Control"):
        # Ladruno (ADR-97 wp/97f, D5): these four are gated behind
        # experimental_integrator now; the baseline itself is untouched
        # (D5 is a parse-time gate, not a behavior change).
        table["cube/vm/%s/Continuum" % meth] = (
            lambda m=meth: cube(lambda t2, m=m: N.mat_vm(t2, "Continuum", method=m,
                                                          experimental=1),
                                N.P_PLASTIC, 6))
    return table[name]()


def child(name, out):
    from _testbed import ops
    rec = {}
    try:
        build, nsteps = _build_one(name)
        build()
    except Exception as exc:
        rec = {"build_error": repr(exc)}
    else:
        codes, sig, eps = [], [], []
        for _ in range(nsteps):
            rc = ops.analyze(1)
            codes.append(int(rc))
            if rc != 0:
                break
            sig.append([float(x) for x in list(ops.eleResponse(1, "stresses"))[0:6]])
            eps.append([float(x) for x in list(ops.eleResponse(1, "strains"))[0:6]])
        rec = {"codes": codes, "stress": sig, "strain": eps}
    with open(out, "w") as f:
        json.dump(rec, f)


def main(out):
    hist = {}
    tmp = out + ".part"
    env = dict(os.environ)
    for name in _decks():
        if os.path.exists(tmp):
            os.remove(tmp)
        p = subprocess.run([sys.executable, os.path.abspath(__file__), "--one", name, tmp],
                           env=env, capture_output=True, text=True)
        if os.path.exists(tmp):
            hist[name] = json.load(open(tmp))
        else:
            hist[name] = {"child_error": (p.stdout + p.stderr)[-400:]}
        v = hist[name]
        print("  %-46s %s" % (name, "steps=%d rc=%s" % (len(v["stress"]), set(v["codes"]))
                              if "stress" in v else list(v)[0]))
        sys.stdout.flush()
    with open(out, "w") as f:
        json.dump(hist, f, indent=1, sort_keys=True)
    n = sum(1 for v in hist.values() if "stress" in v)
    rows = sum(len(v.get("stress", [])) for v in hist.values())
    print("wrote %s: %d decks, %d committed-stress rows" % (out, n, rows))


if __name__ == "__main__":
    if sys.argv[1] == "--one":
        child(sys.argv[2], sys.argv[3])
    else:
        main(sys.argv[1])
