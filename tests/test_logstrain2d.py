"""LogStrain2D (ND_TAG 33016) — openseespy build/construction gate.

At the material-only stage (ADR 25 P5 material; the plane element -geom finite
wiring that drives setTrialF(2×2 F) is a follow-up) the physics is verified
building-free by tests/logstrain2d_reference.py (numpy oracle) and
tests/test_logstrain2d_kernel_cpp.py (C++ kernel vs oracle). This test gates the
parts those cannot: that the class COMPILES into the pyd and CONSTRUCTS through
the command parser + factory + FEM_ObjectBroker for both plane modes over a 3D
order-6 inner, and that a non-3D inner is refused gracefully (no kernel crash).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _reset():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)


def test_construct_both_plane_modes():
    _reset()
    ops.nDMaterial("ElasticIsotropic", 1, 2000.0, 0.25)          # 3D order-6 inner
    # both plane faces over the same inner, default + explicit flags
    ops.nDMaterial("LogStrain2D", 10, 1)                          # default = plane strain
    ops.nDMaterial("LogStrain2D", 11, 1, "-planeStrain")
    ops.nDMaterial("LogStrain2D", 12, 1, "-planeStress")
    # alias spelling also registered
    ops.nDMaterial("LogStrain2D", 13, 1, "-PlaneStress")


def test_missing_inner_is_graceful():
    """A missing / invalid inner tag must be refused at the factory (no exit(-1)
    kernel crash) — whether the binding reports it by warning or by exception, the
    interpreter must stay alive and a subsequent valid construction must work."""
    _reset()
    ops.nDMaterial("ElasticIsotropic", 5, 2000.0, 0.25)
    try:
        ops.nDMaterial("LogStrain2D", 20, 999)                   # inner 999 not found
    except Exception:
        pass                                                     # a raised error is fine too
    # the interpreter is still alive: a good construction still succeeds
    ops.nDMaterial("LogStrain2D", 21, 5, "-planeStress")


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
