"""LadrunoBrick — recorder round-trip + body-force / self-weight validation.

Two of the cheap remaining validation rungs from the plan's "Validation status":

  * RECORDER ROUND-TRIP — the element's setResponse/getResponse stress & strain
    trees must survive the recorder machinery: record via an Element recorder,
    read the file back, and assert it equals the live eleResponse.

  * BODY FORCE / SELF-WEIGHT — the ``-b`` body force and the ``eleLoad
    -selfWeight`` path must integrate to the exact total b*V for EVERY
    formulation. For uri this is the 1-pt centroid integral; for eas it is the
    2x2x2 N-integral added in formEAS (the ``appliedB`` branch), so this is the
    first test that exercises eas body loading.

(The element MASS matrix is validated separately, dynamically, in
test_ladrunoBrick_explicit.py — dt_cr scales like sqrt(rho) and the vibration
period like sqrt(m/k), both gated there.)
"""
import math
import os
import tempfile

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# unit cube (V = 1 exactly), Brick node order: 1-4 = z- face, 5-8 = z+ face
_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]


def _build_cube(formulation, extra=None, E=1000.0, nu=0.3):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    args = ["-formulation", formulation] + (extra or [])
    ops.element("LadrunoBrick", 1, *_CONN, 1, *args)


def _read_last_row(path):
    with open(path) as fh:
        rows = [ln for ln in fh.read().splitlines() if ln.strip()]
    return [float(x) for x in rows[-1].split()]


# --------------------------------------------------------------------------
# 1. recorder round-trip: Element recorder must match the live eleResponse
# --------------------------------------------------------------------------
@pytest.mark.parametrize("formulation", ["std", "bbar", "uri", "eas"])
@pytest.mark.parametrize("resp", ["stresses", "strains"])
def test_element_recorder_roundtrip(formulation, resp):
    """Record `resp` (48 = 8 GP x 6) through an Element recorder; the flushed file
    must equal the live eleResponse for every formulation."""
    _build_cube(formulation)
    # statically-determinate restraints + a uniaxial pull (some real state to record)
    ops.fix(1, 1, 1, 1); ops.fix(2, 0, 1, 1); ops.fix(3, 0, 0, 1); ops.fix(4, 1, 0, 1)
    ops.fix(5, 1, 1, 0); ops.fix(6, 0, 1, 0); ops.fix(8, 1, 0, 0)

    tmp = tempfile.mkdtemp()
    fpath = os.path.join(tmp, f"lb_{formulation}_{resp}.out")
    ops.recorder("Element", "-file", fpath, "-time", "-ele", 1, resp)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, 0.5, 0.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    live = list(ops.eleResponse(1, resp))
    ops.remove("recorders")          # flush the recorder file
    row = _read_last_row(fpath)
    recorded = row[1:]               # drop the leading -time column

    assert len(recorded) == len(live) == 48, f"{resp}: got {len(recorded)} vs live {len(live)}"
    for i, (r, l) in enumerate(zip(recorded, live)):
        tol = 1e-3 + 1e-4 * max(abs(r), abs(l))   # default text-recorder precision
        assert abs(r - l) <= tol, f"{formulation}/{resp}[{i}]: file {r} vs live {l}"


# --------------------------------------------------------------------------
# 2. body force / self-weight integrates to b*V for every formulation
# --------------------------------------------------------------------------
@pytest.mark.parametrize("formulation", ["std", "bbar", "uri", "eas"])
def test_self_weight_total_reaction(formulation):
    """eleLoad -selfWeight (scaling the construction -b) must produce a total base
    reaction equal to the exact body-force resultant b*V. V=1 for the unit cube,
    so the vertical base reaction must equal |bz| for ALL formulations — the
    body-force integral is exact (sum_I N_I = 1) regardless of the rule. This is
    the first gate to exercise eas body loading (formEAS appliedB 2x2x2 loop)."""
    bz = -10.0                       # force per unit volume, downward
    _build_cube(formulation, extra=["-b", 0.0, 0.0, bz])
    for n in (1, 2, 3, 4):           # clamp the z- (base) face fully
        ops.fix(n, 1, 1, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, "-type", "-selfWeight", 0.0, 0.0, 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    ops.reactions()
    Rz = sum(ops.nodeReaction(n, 3) for n in (1, 2, 3, 4))
    assert abs(abs(Rz) - 10.0) <= 1e-6 * 10.0 + 1e-9, (
        f"{formulation}: base vertical reaction {Rz:.8f} should equal |bz|*V=10 "
        f"(body-force integral)"
    )


def test_self_weight_total_is_formulation_independent():
    """The total weight resultant must be identical across all four formulations
    (it depends only on b*V, not on the integration scheme)."""
    bz = -7.0

    def _Rz(formulation):
        _build_cube(formulation, extra=["-b", 0.0, 0.0, bz])
        for n in (1, 2, 3, 4):
            ops.fix(n, 1, 1, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", 1, "-type", "-selfWeight", 0.0, 0.0, 1.0)
        ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0); ops.algorithm("Linear"); ops.analysis("Static")
        assert ops.analyze(1) == 0
        ops.reactions()
        return sum(ops.nodeReaction(n, 3) for n in (1, 2, 3, 4))

    vals = [_Rz(f) for f in ("std", "bbar", "uri", "eas")]
    for v in vals:
        assert abs(abs(v) - 7.0) <= 1e-6 * 7.0 + 1e-9, f"weight {v} != 7"
    assert max(vals) - min(vals) <= 1e-6, f"weight resultant varies across formulations: {vals}"
