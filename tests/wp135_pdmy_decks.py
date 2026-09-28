"""WP-135 -- decks for the PDMY substep-cap fix, shared by the baseline capture
and the pytest battery.

Reuses WP-133's drained plane-strain quad (``wp133_pdmy03_deck``) and adds:

* ``PDMY01`` / ``PDMY02`` versions of the same deck (WP-135 edits all three
  PDMY classes; WP-133's baselines only cover PDMY03), captured from an
  UNMODIFIED build for the byte-identity gate;
* the two-element reproducer of the WP-133 hang (``pair_past_crossing``);
* a one-step "absurd increment" driver (``absurd_step``) that asks a single
  quad for a trial strain far past the substep cap.

Run as a script to capture the byte-identity baseline::

    python -S tests/wp135_pdmy_decks.py <dist/bin> <out.json>

It imports nothing from opensees at module import; callers pass ``ops``.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wp133_pdmy03_deck as D  # noqa: E402

# PDMY01: nd rho G B phi gammaPeak refP d PTAng contrac dilat1 dilat2
#         liquefac1 liquefac2 liquefac3   (then NYS e cs1 cs2 cs3 pa c)
PDMY01 = [2, 2.0, 1.3e5, 2.6e5, 40.0, 0.1, 101.0, 0.5, 26.0,
          0.013, 0.3, 3.0, 5.0, 0.0, 1.0]
# PDMY02: nd rho G B phi gammaPeak refP d PTAng contrac1 contrac3 dilat1 dilat3
#         (then NYS contrac2 dilat2 liquefac1 liquefac2 e cs1 cs2 cs3 pa c)
PDMY02 = [2, 2.0, 1.3e5, 2.6e5, 40.0, 0.1, 101.0, 0.5, 26.0,
          0.013, 0.0, 0.3, 0.0]
# optional tails with c = 0.1 (the WP-133 deck's cohesion; keeps p' > 0)
PDMY01_TAIL = [20, 0.6, 0.9, 0.02, 0.7, 101.0, 0.1]
PDMY02_TAIL = [20, 5.0, 3.0, 1.0, 0.0, 0.6, 0.9, 0.02, 0.7, 101.0, 0.1]

MATS = {
    "PressureDependMultiYield": (PDMY01, PDMY01_TAIL),
    "PressureDependMultiYield02": (PDMY02, PDMY02_TAIL),
    "PressureDependMultiYield03": (D.BASE, D.TAIL),
}


def build(ops, mat, tag, extra=()):
    base, tail = MATS[mat]
    ops.nDMaterial(mat, tag, *base, *tail, *extra)


# PDMY01 on the dense deck is fragile under KrylovNewton: on the UNMODIFIED
# engine it grinds inside analyze() at compression step 16 (a wild iterate,
# |du| ~ 12, then ~1e6 substeps per Gauss point per call -- the WP-135 hang,
# no brake and one element). The byte-identity cases therefore stop at step
# 14, or use a looser sand that ends in an ordinary fast convergence failure.
PDMY01_LOOSE = [2, 2.0, 1.3e5, 2.6e5, 30.0, 0.1, 101.0, 0.5, 26.0,
                0.013, 0.1, 3.0, 5.0, 0.0, 1.0]
PDMY01_GRIND_STEP = 16  # 1-based compression step that ground pre-WP-135


def run_one(ops, mat, nstep=D.NSTEP, base=None):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    if base is None:
        build(ops, mat, 1)
    else:
        ops.nDMaterial(mat, 1, *base, *MATS[mat][1])
    D.set_ele_mat(1, 1)
    D.add_quad(ops, 1, 0, 1, 0.0)
    return D.drive(ops, [(0, 1)], nstep)[1]


CASES = {
    "pdmy01_first14": lambda ops: run_one(ops, "PressureDependMultiYield", 14),
    "pdmy01_loose": lambda ops: run_one(ops, "PressureDependMultiYield",
                                        base=PDMY01_LOOSE),
    "pdmy02": lambda ops: run_one(ops, "PressureDependMultiYield02"),
}


def cs1_mid(ops):
    """WP-133's ``_cs1_mid``: a cs1 that puts the CSL across the path."""
    D.set_ele_mat(1, 1)
    ref = D.run_quad(ops, 1, materials=lambda o: D.build_material(o, 1, tail=D.TAIL))
    k = len(ref["strain"]) // 2
    e = [0.6 + (s[0] + s[1]) * 1.6 for s in ref["strain"]]
    return e[k] + 0.02 * (1000.0 / 101.0) ** 0.7


def pair_past_crossing(ops, cs1, nstep=D.NSTEP):
    """WP-133 G3's two-element pair, driven PAST material 2's crossing.

    Returns (per-element results, rc of the first failed analyze or 0).
    Before WP-135 this stalled inside analyze() on step 20 (the crossing).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    D.build_material(ops, 1, tail=D.TAIL)
    D.build_material(ops, 2, tail=D.TAIL, extra=["-cs1", cs1])
    D.set_ele_mat(1, 1)
    D.set_ele_mat(2, 2)
    D.add_quad(ops, 1, 0, 1, 0.0)
    D.add_quad(ops, 2, 10, 2, 5.0)
    return D.drive(ops, [(0, 1), (10, 2)], nstep)


def absurd_step(ops, mat, ele="quad", dy=-1.0e3):
    """Confine + stage 1 exactly as the deck, then ONE step whose prescribed
    top displacement is ``dy`` on a 1x1 element (|axial strain| = |dy|).

    |dy| = 1e3 asks setSubStrainRate() for ~1e8 substeps per Gauss point
    (octahedral shear / 1e-5), i.e. minutes per material call before WP-135.
    Returns the rc of that analyze(1).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    build(ops, mat, 1)
    D.set_ele_mat(1, 1)
    if ele == "quad":
        D.add_quad(ops, 1, 0, 1, 0.0)
    else:  # SSPquad: a host that SWALLOWS setTrialStrain's return code
        n = [1, 2, 3, 4]
        for k, (x, y) in zip(n, [(0, 0), (1, 0), (1, 1), (0, 1)]):
            ops.node(k, float(x), float(y))
        ops.fix(1, 1, 1)
        ops.fix(2, 0, 1)
        ops.fix(4, 1, 0)
        ops.equalDOF(3, 4, 2)
        ops.element("SSPquad", 1, *n, 1, "PlaneStrain", 1.0)
    old = D.DY
    D.DY = dy
    try:
        out = D.drive(ops, [(0, 1)], 1)
    finally:
        D.DY = old
    return out[1]["n"]  # 1 = the step was accepted, 0 = analyze refused it


def to_hex(res):
    return D.to_hex(res)


def pair_cut(ops, cs1, nstep=D.NSTEP, max_halvings=6):
    """``pair_past_crossing`` with the step-cutting a refusal is meant to
    enable: a failed analyze(1) is retried as 2, 4, ... substeps."""
    orig = ops.analyze
    cuts = [0]

    def analyze(n, *a):
        if n != 1 or not getattr(analyze, "on", False):
            return orig(n, *a)
        rc = orig(1)
        if rc == 0:
            return 0
        # a failed analyze reverts to the last committed state; finish the
        # remaining fraction of the step with halving sub-increments
        remaining, h, rc = 1.0, 0.5, 0
        while remaining > 1e-12:
            h = min(h, remaining)
            oi("LoadControl", D.DY * h)
            if orig(1) == 0:
                remaining -= h
            else:
                cuts[0] += 1
                h *= 0.5
                if h < 0.5 ** max_halvings:
                    rc = -1
                    break
        oi("LoadControl", D.DY)
        return rc

    oi = ops.integrator

    def integrator(*a):
        if a[:2] == ("LoadControl", D.DY):
            analyze.on = True
        return oi(*a)

    ops.analyze, ops.integrator = analyze, integrator
    try:
        out = pair_past_crossing(ops, cs1, nstep)
    finally:
        ops.analyze, ops.integrator = orig, oi
    return out, cuts[0]


def _scenario(ops, name):
    """One named scenario, for the subprocess (wall-clock guarded) gates."""
    if name == "pair":
        out = pair_past_crossing(ops, cs1_mid(ops))
        return {"n1": out[1]["n"], "n2": out[2]["n"]}
    if name == "pair_cut":
        out, cuts = pair_cut(ops, cs1_mid(ops))
        return {"n1": out[1]["n"], "n2": out[2]["n"], "cuts": cuts}
    if name == "pdmy01_full":
        return {"n": run_one(ops, "PressureDependMultiYield")["n"]}
    if name.startswith("absurd:"):
        _, mat, ele = name.split(":")
        return {"n": absurd_step(ops, mat, ele)}
    raise ValueError(name)


def _bind(dist):
    os.add_dll_directory(dist)
    os.environ["PATH"] = dist + os.pathsep + os.environ.get("PATH", "")
    sys.path.insert(0, dist)
    import opensees as ops
    assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(dist), ops.__file__
    return ops


def main(dist, out):
    ops = _bind(dist)
    data = {name: to_hex(fn(ops)) for name, fn in CASES.items()}
    with open(out, "w") as f:
        json.dump(data, f, indent=0)
    for k, v in data.items():
        print(k, v["n"], v["stress"][-1] if v["n"] else None)


if __name__ == "__main__":
    if sys.argv[2] == "run":  # <dist> run <scenario> -> one RESULT json line
        res = _scenario(_bind(sys.argv[1]), sys.argv[3])
        print("RESULT " + json.dumps(res), flush=True)
    else:
        main(sys.argv[1], sys.argv[2])
