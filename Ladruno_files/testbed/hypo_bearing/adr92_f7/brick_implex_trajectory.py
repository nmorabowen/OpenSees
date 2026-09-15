"""WP-99 (F7) -- the LadrunoBrick IMPL-EX trajectory recorder.

Writes one CSV row per committed step of a free-DOF LadrunoBrick drained
triaxial under `-implex` with an ADEQUATE `-maxSubsteps` cap (20000), i.e. a
deck on which the commit-time companion cap is NEVER hit and WP-99's latch is
therefore never armed.  Run it against the PRE-change and POST-change binaries
and `fc` / `cmp` the two files: they must be byte-identical.  That is the
"LadrunoBrick bit-identical" claim, diffed rather than asserted.

Usage (from the worktree root):

    PYTHONPATH=<worktree>/dist/bin python3.12 \
        Ladruno_files/testbed/hypo_bearing/adr92_f7/brick_implex_trajectory.py out.csv

The deck is `tests/test_ladruno_sanisand_implex.py::_build_free_dof_triaxial`
in shape (the same one the refusal tests use) with the CP1/Gorini parameter
set, kept here rather than imported so the recorder runs standalone.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "tests"))

import opensees as ops   # noqa: E402  (PYTHONPATH must point at dist/bin)

# CP1 / Gorini parameter set -- tests/test_ladruno_sanisand.py::_PARAMS
_PARAMS = [
    264.32, 0.3129, 0.6944, 1.33090, 0.71, 0.027, 0.83, 0.45, 101.0,
    0.005, 1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0,
]
_P_ATM = _PARAMS[8]
_XY = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]

_P0 = 100.0            # kPa
_N_CONF = 5
_N_DEV = 12
_DQ = 6.0              # kPa per deviatoric step
_TOL_REL = 1.0e-3
_MAXITER = 60
_CAP_ADEQUATE = 20000


def _build(tag):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial("LadrunoSANISAND", tag, *_PARAMS,
                   1, 2, 1, 1.0e-7, 1.0e-7,
                   "-Presidual", 0.0, "-Pmin", 1.0e-4 * _P_ATM,
                   "-implex", "-maxSubsteps", _CAP_ADEQUATE)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, tag,
                "-geom", "linear", "-formulation", "bbar")
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    q = _P0 / 4.0
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.load(n, -q, 0.0, 0.0)
            if y == 1.:
                ops.load(n, 0.0, -q, 0.0)
            if k == 1:
                ops.load(n, 0.0, 0.0, -q)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormUnbalance", _TOL_REL * _P0, _MAXITER, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / _N_CONF)
    ops.analysis("Static")


def _row(step):
    sig = list(ops.eleResponse(1, "material", 1, "stress"))
    eps = list(ops.eleResponse(1, "material", 1, "strain"))
    det = list(ops.eleResponse(1, "material", 1, "implexDetail"))
    ref = list(ops.eleResponse(1, "material", 1, "implexRefusals"))
    vals = sig + eps + det[:6] + ref[:4]
    return ",".join([str(step)] + ["%.17g" % v for v in vals])


def main(out_path):
    tag = 9901
    _build(tag)
    rows = []

    ops.updateMaterialStage("-material", tag, "-stage", 0)
    for step in range(_N_CONF):
        assert ops.analyze(1) == 0, "confinement step %d failed" % (step + 1)
        rows.append(_row("conf%d" % (step + 1)))
    ops.loadConst("-time", 0.0)
    ops.updateMaterialStage("-material", tag, "-stage", 1)

    dq = _DQ / 4.0
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -dq)
    ops.integrator("LoadControl", 1.0)
    for step in range(_N_DEV):
        assert ops.analyze(1) == 0, "deviatoric step %d failed" % (step + 1)
        rows.append(_row("dev%d" % (step + 1)))

    header = ("step,"
              + ",".join("sig%d" % i for i in range(6)) + ","
              + ",".join("eps%d" % i for i in range(6)) + ","
              + "err_total,err_dev,err_vol,clampFired,clampCount,f,"
              + "ref_total,ref_signChange,ref_control,ref_companion")
    with open(out_path, "w", newline="\n") as fh:
        fh.write(header + "\n")
        for r in rows:
            fh.write(r + "\n")
    print("wrote %s (%d rows)" % (out_path, len(rows)))
    print("build:", ops.ladrunoBuild().strip().splitlines()[0])


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "brick_implex_trajectory.csv")
