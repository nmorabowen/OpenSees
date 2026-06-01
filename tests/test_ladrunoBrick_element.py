"""LadrunoBrick (Ladruno unified 8-node hexahedron, classTag 33002) — Zone-A
element battery.

v1 ships the small-strain ``std`` + ``bbar`` formulations. The headline gate is
**bit-for-bit overlap with the upstream elements**:

  * ``-formulation std``  must reproduce ``stdBrick``   (Ed Love's Brick)
  * ``-formulation bbar`` must reproduce ``bbarBrick``  (Ed Love's BbarBrick)

to ~1e-9 on a distorted hex under mixed loading. That equivalence is our
correctness anchor (only possible because v1 is geometrically linear). On top of
that:

  * rank / zero-energy: a free element has exactly 6 rigid-body modes (18 with
    energy) for both std and bbar — no spurious mechanisms;
  * volumetric-locking relief: near-incompressible (nu=0.4999), bbar is strictly
    more flexible than std (std locks);
  * the reserved formulations (uri/eas) are refused at the factory in v1.

Plan: Ladruno_implementation/09_ladruno_brick.md.
"""
import math

import pytest

from _testbed import ops
from _testbed.fem_checks import zero_energy_mode_count, n_rigid_body_modes

pytestmark = [pytest.mark.zone_a]

# A deliberately distorted (but positive-Jacobian) hex. Standard Brick node
# order: nodes 1-4 are the z- face (CCW), nodes 5-8 the z+ face (CCW).
_NODES = {
    1: (0.00, 0.00, 0.00),
    2: (1.00, 0.10, 0.00),
    3: (1.10, 1.00, 0.10),
    4: (0.05, 0.95, 0.00),
    5: (0.00, 0.05, 1.00),
    6: (1.00, 0.00, 1.05),
    7: (1.05, 1.00, 1.10),
    8: (0.00, 1.00, 0.95),
}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]
_BASE = [1, 2, 3, 4]          # fully fixed
_TOP = [5, 6, 7, 8]           # loaded
_LOADS = {                    # mixed so the full 24x24 K is exercised
    5: (10.0, 0.0, 0.0),
    6: (0.0, 8.0, 0.0),
    7: (0.0, 0.0, -12.0),
    8: (5.0, 5.0, 5.0),
}


def _build_common(E, nu):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)


def _solve_hex(ele_name, extra_args, E=1000.0, nu=0.3):
    """Build a single hex of the given element, static-solve under _LOADS, and
    return (disps[24], stresses[48], forces[24])."""
    _build_common(E, nu)
    for n in _BASE:
        ops.fix(n, 1, 1, 1)

    ops.element(ele_name, 1, *_CONN, 1, *extra_args)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (fx, fy, fz) in _LOADS.items():
        ops.load(n, fx, fy, fz)

    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    disps = []
    for n in _CONN:
        disps.extend(ops.nodeDisp(n))
    # Read "forces" BEFORE "stresses": upstream bbarBrick has no update() and
    # sets the material trial strain only inside formResidAndTangent, so a bare
    # "stresses" read reflects the predictor (u=0) state. getResistingForce
    # ("forces") triggers that strain evaluation at the committed displacement,
    # after which "stresses" reflects the solved state for lazy-strain elements
    # too. (LadrunoBrick/Brick implement update(), so order is immaterial for
    # them.) See LEDGER_quirks: bbarBrick lazy-strain readback.
    forces = list(ops.eleResponse(1, "forces"))
    stresses = list(ops.eleResponse(1, "stresses"))
    return disps, stresses, forces


def _assert_close(a, b, label, rtol=1e-9, atol=1e-9):
    assert len(a) == len(b), f"{label}: length mismatch {len(a)} vs {len(b)}"
    for i, (x, y) in enumerate(zip(a, b)):
        tol = atol + rtol * max(abs(x), abs(y))
        assert abs(x - y) <= tol, (
            f"{label}[{i}]: {x!r} vs {y!r} (|d|={abs(x-y):.3e} > {tol:.3e})"
        )


# --------------------------------------------------------------------------
# headline regression gates
# --------------------------------------------------------------------------
def test_std_matches_stdBrick():
    """-formulation std reproduces upstream stdBrick to ~1e-9."""
    ref = _solve_hex("stdBrick", [])
    ours = _solve_hex("LadrunoBrick", ["-formulation", "std"])
    _assert_close(ours[0], ref[0], "disp")
    _assert_close(ours[1], ref[1], "stress")
    _assert_close(ours[2], ref[2], "force")


def test_bbar_matches_bbarBrick():
    """-formulation bbar reproduces upstream bbarBrick to ~1e-9."""
    ref = _solve_hex("bbarBrick", [])
    ours = _solve_hex("LadrunoBrick", ["-formulation", "bbar"])
    _assert_close(ours[0], ref[0], "disp")
    _assert_close(ours[1], ref[1], "stress")
    _assert_close(ours[2], ref[2], "force")


def test_default_formulation_is_std():
    """No -formulation flag defaults to std."""
    ref = _solve_hex("LadrunoBrick", ["-formulation", "std"])
    bare = _solve_hex("LadrunoBrick", [])
    _assert_close(bare[0], ref[0], "disp")


# --------------------------------------------------------------------------
# rank / zero-energy (no spurious mechanisms for std/bbar)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("formulation", ["std", "bbar"])
def test_free_element_has_six_rigid_body_modes(formulation):
    _build_common(E=1000.0, nu=0.3)
    ops.element("LadrunoBrick", 1, *_CONN, 1, "-formulation", formulation)
    count, eigs = zero_energy_mode_count(_CONN, ndf=3)
    assert count == n_rigid_body_modes(3), (
        f"{formulation}: {count} zero-energy modes (want 6); eigs={eigs}"
    )


# --------------------------------------------------------------------------
# volumetric-locking relief: bbar is strictly more flexible near-incompressible
# --------------------------------------------------------------------------
def test_bbar_relieves_volumetric_locking():
    nu = 0.4999
    std = _solve_hex("LadrunoBrick", ["-formulation", "std"], nu=nu)
    bbar = _solve_hex("LadrunoBrick", ["-formulation", "bbar"], nu=nu)
    std_mag = math.sqrt(sum(u * u for u in std[0]))
    bbar_mag = math.sqrt(sum(u * u for u in bbar[0]))
    # std locks (over-stiff) -> smaller displacement magnitude than bbar.
    assert bbar_mag > std_mag, (
        f"bbar (|u|={bbar_mag:.3e}) should be more flexible than std "
        f"(|u|={std_mag:.3e}) at nu={nu}"
    )


# --------------------------------------------------------------------------
# reserved formulations refused in v1
# --------------------------------------------------------------------------
@pytest.mark.parametrize("formulation", ["uri", "eas"])
def test_reserved_formulations_refused(formulation):
    _build_common(E=1000.0, nu=0.3)
    try:
        ops.element("LadrunoBrick", 1, *_CONN, 1, "-formulation", formulation)
    except Exception:
        pass  # raising is an acceptable refusal
    tags = ops.getEleTags() or []
    if isinstance(tags, int):
        tags = [tags]
    assert 1 not in tags, f"-formulation {formulation} should be refused in v1"
