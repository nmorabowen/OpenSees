"""WP-232l T1: route E equals LadrunoUP (tests/data/mpm/twophase/t1_fe.json).

Linear basis, particles at the Gauss points with V_p = the Gauss weight x detJ, `kinematics
small`, from the common state of case t1_ref_fe_identity (sigma' = 0, F = I, prescribed particle
v, a, p, pd, pdd; the FE nodal state is its P5 projection). `mpm debug assemble` gives the
Newmark first-iteration tangent A and unbalance b: A within rel 1e-12 (Frobenius) and b within
rel 1e-12 (2-norm) of LadrunoUP's; the first step's nodal u and p within rel 1e-10. Twelve legs:
Q4 one cell and 2 x 2, H8 one cell; ElasticIsotropic and DruckerPrager; stab off and auto 0.25.
The q4_1 elastic legs are also the T11 Zone-A oracle: skipping Q or the p projection turns them red.
"""
from __future__ import annotations

import os
import sys
import tempfile

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_twophase_harness as H  # noqa: E402

pytestmark = pytest.mark.zone_a

LEGS = ["q4_1_elastic_off", "q4_1_elastic_auto", "q4_1_dp_off", "q4_1_dp_auto",
        "q4_2x2_elastic_off", "q4_2x2_elastic_auto", "q4_2x2_dp_off", "q4_2x2_dp_auto",
        "h8_1_elastic_off", "h8_1_elastic_auto", "h8_1_dp_off", "h8_1_dp_auto"]


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


@pytest.mark.parametrize("leg", LEGS)
def test_t1_fe_identity(leg):
    ops = _ops()
    rec = H.load(H.T1_REC)
    with tempfile.TemporaryDirectory() as tmp:
        row = H.run_t1(ops, rec, tmp, legs=[leg])[0]
    assert H.t1_pass(row), row
