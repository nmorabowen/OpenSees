"""A refused system, numberer, constraints, test or integrator command is an
error and leaves the previous object in force.

Each case defines a complete set of analysis objects, issues one command that
creates nothing (an unknown integrator name, or a known name whose factory
rejects its options), and then rebuilds the analysis from the stored objects.

Model: one truss, E = 1, A = 1, L = 1, load P = 1 at node 2, so k = EA/L = 1.
The previous objects are LoadControl 0.1 and FixedNumIter 3. One step must
give t = 0.1, u2 = 0.1 * P / k = 0.1 and testIter = 3. A stored object that was
nulled is replaced by a default when the analysis is rebuilt: the analysis
command prints "no LinearSOE specified" (or the numberer, handler, integrator
equivalent), the integrator default LoadControl 1.0 gives t = 1.0, and the
test default NormUnbalance 1e-6 25 converges this linear model in 1 iteration.

Before the fix the refused commands returned 0, so OpenSeesPy raised nothing,
and system, numberer, constraints and test installed a null object.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


DT = 0.1
DEFAULT_WARNINGS = (
    "no LinearSOE specified",
    "no Numberer specified",
    "no ConstraintHandler yet specified",
    "no Integrator specified",
)


def _model():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("FixedNumIter", 3)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", DT)


# (command, arguments for which no object is created)
REFUSED = [
    pytest.param("system", ("PFEM", "-noSuchOption"), id="system-PFEM-unknown-option"),
    pytest.param("numberer", ("ParallelRCM",), id="numberer-ParallelRCM-serial"),
    pytest.param("constraints", ("Penalty",), id="constraints-Penalty-no-alpha"),
    pytest.param("test", ("NormUnbalance",), id="test-NormUnbalance-no-tol"),
    pytest.param("integrator", ("LoadControl",), id="integrator-LoadControl-no-dlam"),
    pytest.param("integrator", ("NoSuchIntegrator", 1.0), id="integrator-unknown-name"),
]


@pytest.mark.parametrize("cmd,refused", REFUSED)
@pytest.mark.parametrize("analysis_first", [False, True],
                         ids=["no-analysis", "with-analysis"])
def test_refused_object_raises_and_keeps_previous(cmd, refused, analysis_first, capfd):
    try:
        _model()
        if analysis_first:
            ops.analysis("Static")
        capfd.readouterr()
        with pytest.raises(Exception):
            getattr(ops, cmd)(*refused)
        assert "left unchanged" in capfd.readouterr().err

        ops.analysis("Static")
        err = capfd.readouterr().err
        for warning in DEFAULT_WARNINGS:
            assert warning not in err, err
        assert ops.analyze(1) == 0
        assert ops.getTime() == pytest.approx(DT, rel=1e-12)
        assert ops.nodeDisp(2, 1) == pytest.approx(DT, rel=1e-12)
        assert ops.testIter() == 3
    finally:
        ops.wipe()
