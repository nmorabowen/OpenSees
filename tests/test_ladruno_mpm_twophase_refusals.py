"""WP-232l T8 (part): two-phase refusals and options, mode A (interface P7, §7a Q4).

- `mpm massCutoff` leaves route S's `mpm pressureCutoff` unchanged (ruling c84f71fb3), and
  `mpm info twoPhase` reports it.
- Loud refusals at `mpm analyze`: a symmetric-storage system, mode quasiStatic, HHT, an
  initial-stress offset on a two-phase particle, `mpm antiLocking fbar`, `-dynSeepage on`.
- Mode F refuses a two-phase body until route S lands (nothing is ignored silently).
"""
from __future__ import annotations

import pytest

pytestmark = pytest.mark.zone_a


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


def _as_list(v):
    return [float(x) for x in v] if isinstance(v, (list, tuple)) else [float(v)]


def _deck(ops, mode="implicit", stress=False, extra=()):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 1.0, "-cells", 2, 2, "-basis", "gimp")
    ops.mpm("mode", mode)
    ops.mpm("body", 1, "box", 0.0, 0.0, 2.0, 2.0, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    if stress:
        ops.mpm("stress", 1, "uniform", -1.0, -1.0, 0.0)
    ops.mpm("twoPhase", 1, "-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-perm", 0.01, 0.01)
    ops.mpm("fix", "ymin", "y")
    ops.mpm("pfix", "ymax", 0.0)
    for a in extra:
        ops.mpm(*a)


def _analysis(ops, system="UmfPack", integrator=("Newmark", 0.6, 0.3025), static=False):
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system(system)
    ops.algorithm("Linear")
    if static:
        ops.integrator("LoadControl", 0.1)
        ops.analysis("Static")
    else:
        ops.integrator(*integrator)
        ops.analysis("Transient")


def _refused(ops):
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 0.01)


def test_pressure_cutoff_is_its_own_knob():
    ops = _ops()
    _deck(ops)
    assert _as_list(ops.mpm("info", "twoPhase", 1))[-1] == 5e-2
    ops.mpm("massCutoff", 1e-10)
    assert _as_list(ops.mpm("info", "twoPhase", 1))[-1] == 5e-2
    ops.mpm("pressureCutoff", 0.1)
    assert _as_list(ops.mpm("info", "twoPhase", 1))[-1] == 0.1
    ops.mpm("massCutoff", 0.3)
    assert _as_list(ops.mpm("info", "twoPhase", 1))[-1] == 0.1


def test_admitted_deck_runs():
    ops = _ops()
    _deck(ops)
    _analysis(ops)
    ops.mpm("analyze", 1, 0.01)
    p = _as_list(ops.mpm("particleData", "p"))
    assert len(p) == 16 and all(x == x for x in p)


@pytest.mark.parametrize("why", ["symmetric", "static", "hht", "sig0", "fbar"])
def test_refusals(why):
    ops = _ops()
    _deck(ops, mode="quasiStatic" if why == "static" else "implicit", stress=(why == "sig0"),
          extra=([("antiLocking", "fbar")] if why == "fbar" else []))
    if why == "symmetric":
        _analysis(ops, system="ProfileSPD")
    elif why == "static":
        _analysis(ops, static=True)
    elif why == "hht":
        _analysis(ops, integrator=("HHT", 0.9))
    else:
        _analysis(ops)
    _refused(ops)


def test_dynseepage_on_refused():
    ops = _ops()
    _deck(ops)
    with pytest.raises(Exception):
        ops.mpm("twoPhase", 1, "-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-perm", 0.01, 0.01,
                "-dynSeepage", "on")


def test_mode_f_refuses_until_route_s():
    ops = _ops()
    _deck(ops, mode="explicit")
    _refused(ops)
