"""Self-test for quirk-lint rule `ground-sign` (alias L6, WP-117) in ci/check_quirk_patterns.py.

One file per rule (WP-162), so two PRs adding rules never edit the same test file.
Shared helpers live in ci/_quirk_testkit.py.
Run: pytest -q ci/test_quirk_ground_sign.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402


RF_SUB_Q = """
const Vector &Elem::getResistingForce()
{
    P_return.Zero();
    P_return.addVector(1.0, Q, -1.0);
    return P_return;
}
"""
BEZIER_WRONG = RF_SUB_Q + """
int Elem::addInertiaLoadToUnbalance(const Vector &accel)
{
    const Matrix &M = this->getMass();
    static Vector a(6);
    // Q += M x a
    Q.addMatrixVector(1.0, M, a, 1.0);
    return 0;
}
"""


def _run(tmp_path, body):
    root = kit.tree(tmp_path, {"SRC/element/Elem.cpp": body})       # unstamped: L6 scans vanilla too
    return cq.check_ground_sign(root, kit.rel(root), set())


def test_flags_the_bezier_shape(tmp_path):
    out = _run(tmp_path, BEZIER_WRONG)
    assert len(out) == 1 and "wrong sign" in out[0] and "+M*R*a_g" in out[0], out


def test_passes_the_wp117_fix(tmp_path):
    assert _run(tmp_path, BEZIER_WRONG.replace("Q.addMatrixVector(1.0, M, a, 1.0);",
                                              "Q.addMatrixVector(1.0, M, a, -1.0);")) == []


def test_reads_the_other_accumulation_forms(tmp_path):
    # FourNodeQuad: `for (...) Q(i) += -K(i,i)*ra[i];` (for-header semicolons must not split it)
    quad = RF_SUB_Q + """
int Elem::addInertiaLoadToUnbalance(const Vector &accel)
{
  for (int i = 0; i < 8; i++)
    Q(i) += -K(i,i)*ra[i];
  return 0;
}
"""
    assert _run(tmp_path / "a", quad) == []
    assert len(_run(tmp_path / "b", quad.replace("Q(i) += -K(i,i)*ra[i];", "Q(i) += K(i,i)*ra[i];"))) == 1
    # ElasticBeam: `Q(0) -= m * Raccel1(0);`
    beam = RF_SUB_Q + """
int Elem::addInertiaLoadToUnbalance(const Vector &accel)
{
  Q(0) -= m * Raccel1(0);
  Q(1) -= m * Raccel1(1);
  return 0;
}
"""
    assert _run(tmp_path / "c", beam) == []
    # LadrunoBrick: `load->addMatrixVector(1.0, M, resid, -1.0);` with `resid -= *load;`
    brick = """
const Vector &Elem::getResistingForce(void)
{
  formResidAndTangent(0);
  if (load != 0) resid -= *load;
  return resid;
}
int Elem::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (load == 0) load = new Vector(24);
  load->addMatrixVector(1.0, M, resid, -1.0);
  return 0;
}
"""
    assert _run(tmp_path / "d", brick) == []
    assert len(_run(tmp_path / "e", brick.replace("resid, -1.0)", "resid, 1.0)"))) == 1


def test_added_vector_needs_a_positive_accumulation(tmp_path):
    # the other consistent convention: accumulate +M*a and ADD it to the residual
    added = BEZIER_WRONG.replace("P_return.addVector(1.0, Q, -1.0);", "P_return.addVector(1.0, Q, 1.0);")
    assert _run(tmp_path, added) == []


def test_stays_silent_when_the_sign_is_unreadable(tmp_path):
    # a variable factor cannot be judged: skip, never guess
    assert _run(tmp_path, BEZIER_WRONG.replace("a, 1.0);", "a, fact);")) == []


def test_waiver(tmp_path):
    waived = BEZIER_WRONG.replace("    Q.addMatrixVector(1.0, M, a, 1.0);",
                                  "    // ladruno-lint: sign-ok Q here is ADDED by a custom integrator hook\n"
                                  "    Q.addMatrixVector(1.0, M, a, 1.0);")
    assert _run(tmp_path, waived) == []
