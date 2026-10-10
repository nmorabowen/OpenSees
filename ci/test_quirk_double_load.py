"""Self-test for quirk-lint rule `double-load` (alias L5, WP-119) in ci/check_quirk_patterns.py.

One file per rule (WP-162), so two PRs adding rules never edit the same test file.
Shared helpers live in ci/_quirk_testkit.py.
Run: pytest -q ci/test_quirk_double_load.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402


# the vanilla ElasticBeam2d shape (pre WP-119): RF subtracts Q, IncInertia calls RF and subtracts Q again
RF_SUBTRACTS_Q = """
const Vector &
Beam::getResistingForce()
{
  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);
  if (rho != 0)
    P.addVector(1.0, Q, -1.0);
  return P;
}
"""
DOUBLE_Q = RF_SUBTRACTS_Q + """
const Vector &
Beam::getResistingForceIncInertia()
{
  P = this->getResistingForce();
  // subtract external load P = P - Q
  P.addVector(1.0, Q, -1.0);
  if (rho == 0.0)
    return P;
  P(1) += m * accel1(1);
  return P;
}
"""


def _run(tmp_path, body, where="SRC/element/Beam.cpp"):
    root = kit.tree(tmp_path, {where: body})
    return cq.check_double_load(root, kit.rel(root), set())


def test_flags_the_elasticbeam2d_shape_in_a_vanilla_file(tmp_path):
    out = _run(tmp_path, DOUBLE_Q)                          # no LADRUNO stamp: vanilla is in scope for L5
    assert len(out) == 1 and "subtracts 'Q' again" in out[0], out


def test_passes_the_wp119_fix(tmp_path):
    fixed = DOUBLE_Q.replace("  P.addVector(1.0, Q, -1.0);\n  if (rho == 0.0)",
                             "  // P.addVector(1.0, Q, -1.0);\n  if (rho == 0.0)")
    assert _run(tmp_path, fixed) == []


def test_flags_the_timoshenko_shape(tmp_path):
    timo = """
const Vector& Beam::getResistingForce()
{
    theVector.addMatrixTransposeVector(0.0, Tgl, ql, 1.0);
    if (rho != 0.0)
      theVector.addVector(1.0, theLoad, -1.0);
    return theVector;
}
const Vector& Beam::getResistingForceIncInertia()
{
    theVector = this->getResistingForce();
    theVector.addVector(1.0, theLoad, -1.0);
    return theVector;
}
"""
    assert len(_run(tmp_path, timo)) == 1


def test_passes_the_single_subtraction_shapes(tmp_path):
    # FourNodeQuad shape: RF subtracts Q, IncInertia calls RF and only adds inertia
    quad = RF_SUBTRACTS_Q + """
const Vector &
Beam::getResistingForceIncInertia()
{
  P = this->getResistingForce();
  P.addMatrixVector(1.0, M, a, 1.0);
  return P;
}
"""
    assert _run(tmp_path / "a", quad) == []
    # LadrunoBrick shape: IncInertia rebuilds the residual without calling getResistingForce()
    brick = RF_SUBTRACTS_Q + """
const Vector &
Beam::getResistingForceIncInertia()
{
  formResidAndTangent(0);
  res = resid;
  if (load != 0) res -= *load;
  return res;
}
"""
    assert _run(tmp_path / "b", brick) == []
    # a different vector subtracted in IncInertia (e.g. an initial-force offset) is not the load
    other = DOUBLE_Q.replace("  P.addVector(1.0, Q, -1.0);\n  if (rho == 0.0)",
                             "  P.addVector(1.0, P0, -1.0);\n  if (rho == 0.0)")
    assert _run(tmp_path / "c", other) == []


def test_waiver(tmp_path):
    waived = DOUBLE_Q.replace("  P.addVector(1.0, Q, -1.0);\n  if (rho == 0.0)",
                              "  // ladruno-lint: double-ok Q here is a distinct follower-load buffer\n"
                              "  P.addVector(1.0, Q, -1.0);\n  if (rho == 0.0)")
    assert _run(tmp_path, waived) == []
