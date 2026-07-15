"""Regression: element-level `rho` must survive a database/parallel restart.

FourNodeQuad and its siblings (Tri31, Nine/Eight/SixNode, FourNodeQuad3d) carry
an element-level density `rho` that overrides the material density in the mass
matrix via `if (rho == 0) use material->getRho() else use rho`. Upstream never
serialized `rho` in sendSelf/recvSelf, and the broker (no-arg) constructor left
it UNINITIALIZED — so after `database`/`restore` (or an OpenSeesMP send) the
element's `rho` was garbage heap memory. Garbage != 0, so it hijacked the mass
matrix: non-deterministic, indefinite M (negative generalized eigenvalues),
diverging transient restarts. See LEDGER_quirks.md "element rho un-serialized".

The committed nodal disp round-trips fine (that is why this was invisible to a
plain `database_roundtrip`); the corruption only shows once the CONTINUED
transient uses the mass matrix. So this test continues the integration across the
restart and compares to the uninterrupted run.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

DT = 0.01
NY = 4  # node rows (2 wide)


def _nid(col, row):
    return row * 2 + col + 1


def _build(elem):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for row in range(NY):
        for col in (0, 1):
            ops.node(_nid(col, row), float(col), float(row))
    ops.fix(_nid(0, 0), 1, 1)
    ops.fix(_nid(1, 0), 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)  # rho=2 lives on the material
    et = 1
    if elem == "quad":
        for r in range(NY - 1):
            ops.element("quad", et, _nid(0, r), _nid(1, r), _nid(1, r + 1),
                        _nid(0, r + 1), 1.0, "PlaneStrain", 1)
            et += 1
    elif elem == "tri31":
        for r in range(NY - 1):
            ops.element("tri31", et, _nid(0, r), _nid(1, r), _nid(1, r + 1),
                        1.0, "PlaneStrain", 1); et += 1
            ops.element("tri31", et, _nid(0, r), _nid(1, r + 1), _nid(0, r + 1),
                        1.0, "PlaneStrain", 1); et += 1
    else:
        raise ValueError(elem)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(_nid(0, NY - 1), 0.0, 1.0)  # excitation so the mass term matters


def _setup():
    ops.constraints("Penalty", 1e12, 1e12)
    ops.numberer("Plain")
    ops.system("ProfileSPD")
    ops.test("NormDispIncr", 1e-10, 20, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")


def _top():
    return ops.nodeDisp(_nid(0, NY - 1), 2)


@pytest.mark.parametrize("elem", ["quad", "tri31"])
def test_element_rho_survives_db_restart(elem, tmp_path):
    # uninterrupted reference: 6 continuous steps
    _build(elem)
    _setup()
    ref = []
    for _ in range(6):
        assert ops.analyze(1, DT) == 0
        ref.append(_top())
    assert max(abs(v) for v in ref) > 1e-6, "reference state ~zero; test would be vacuous"

    # restart: 3 steps, save -> wipe -> restore, then 3 more steps
    _build(elem)
    _setup()
    for _ in range(3):
        assert ops.analyze(1, DT) == 0
    dbpath = str(tmp_path / "rho_rt")
    try:
        ops.database("File", dbpath)
    except Exception as exc:  # noqa: BLE001 - build without FE_Datastore
        pytest.skip(f"database() unsupported in this build: {exc}")
    saved = ops.save(1)
    if saved is not None and saved < 0:
        pytest.skip("database save returned failure on this build")
    ops.wipe()
    ops.database("File", dbpath)
    ops.restore(1)
    _setup()
    got = []
    for _ in range(3):
        assert ops.analyze(1, DT) == 0
        got.append(_top())
    ops.wipe()  # release FE_Datastore handles (Windows) before tmp cleanup

    # steps 4-6 of the restart must reproduce the uninterrupted run; a corrupt
    # (un-restored) element rho gives an indefinite mass matrix and gross drift.
    for i, (a, b) in enumerate(zip(got, ref[3:])):
        assert abs(a - b) <= 1e-6 * abs(b) + 1e-12, (
            f"{elem}: DB restart step {i + 4} diverged {b} -> {a}; element rho "
            "not serialized/restored -> corrupt mass matrix"
        )
