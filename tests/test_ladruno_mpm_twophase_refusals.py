"""WP-232l T8 (part): two-phase refusals and options, mode A (interface P7, §7a Q4).

- `mpm massCutoff` leaves route S's `mpm pressureCutoff` unchanged (ruling c84f71fb3), and
  `mpm info twoPhase` reports it.
- Loud refusals at `mpm analyze`: a symmetric-storage system, mode quasiStatic, HHT, an
  initial-stress offset on a two-phase particle, `mpm antiLocking fbar`, `-dynSeepage on`.
- Mode F refuses a two-phase body until route S lands (nothing is ignored silently).
"""
from __future__ import annotations

import math

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


def test_generalized_alpha_refused():
    ops = _ops()
    _deck(ops)
    _analysis(ops, integrator=("GeneralizedAlpha", 0.9, 0.8))
    _refused(ops)


@pytest.mark.parametrize("args", [
    ["-poro", 0.4, "-rhoF", 1.0, "-perm", 0.01, 0.01],                         # no -Kf
    ["-Kf", 4444.0, "-poro", 0.6, "-rhoF", 1.0, "-perm", 0.01, 0.01, "-alpha", 0.5],   # n > alpha
    ["-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-perm", 0.01, 0.01, "-bogus", 1.0],   # unknown flag
    ["-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-permH", 0.01, 0.01],          # -permH without -gammaW
])
def test_twophase_parser_refusals(args):
    ops = _ops()
    _deck(ops)
    with pytest.raises(Exception):
        ops.mpm("twoPhase", 2, *args)


def test_porous_material_refused_at_particles():
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    ops.nDMaterial("FluidSolidPorous", 2, 2, 1, 2.2e6, 101.0)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 1.0, "-cells", 2, 2)
    with pytest.raises(Exception):
        ops.mpm("body", 1, "box", 0.0, 0.0, 2.0, 2.0, "-ppc", 2, "-mat", 2, "-rho", 2.0)


def test_bspline_pfix_nonzero_and_traction_refused():
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 1.0, "-cells", 2, 2, "-basis", "bspline2")
    ops.mpm("pfix", "ymax", 0.0)                     # drained: admitted
    with pytest.raises(Exception):
        ops.mpm("pfix", "xmax", 3.0)
    with pytest.raises(Exception):
        ops.mpm("traction", "t", "ymax", 0.0, 2.0, 1.0, 2.0, 0.0, -1.0)


def test_pack_below_5_refused_with_two_phase():
    ops = _ops()
    _deck(ops)
    with pytest.raises(Exception):
        ops.mpm("debug", "packVersion", 4)


def test_info_reports_route_constants_and_stab():
    ops = _ops()
    _deck(ops)
    info = _as_list(ops.mpm("info", "twoPhase", 1))
    # [route, Kf, n, rhoF, k_x, k_y, alpha, Ks, stabMode, a0, alpha_st(h), pressureCutoff]
    assert info[:10] == [1.0, 4444.0, 0.4, 1.0, 0.01, 0.01, 1.0, -1.0, 1.0, 0.25]
    G = 1.0e4 / 2.4
    lam = 1.0e4 * 0.2 / (1.2 * 0.6)
    assert info[10] == pytest.approx(0.25 / (lam + 2 * G / 3 + 4 * G / 3), rel=1e-12)


@pytest.mark.xfail(strict=True, reason=(
    "WP-232l T8 miss on fd81cbdd3 (reported, not fixed in this binary): P4's literal stray-p rule "
    "(p-active iff a two-phase particle has S^p > 0 there) admits a node that a two-phase particle "
    "touches with a GIMP tail of ~1e-8 after one step of motion; that p row is near-singular and "
    "p reaches -6.7e3 at step 2. Proposed fix: mode A's p-activity on the two-phase nodal mass "
    "against mode A's cutoff, S^p renormalised over p-active nodes (route S's ruling 3d3699a50)."))
def test_stray_p_is_fixed_where_no_two_phase_particle_reaches():
    """U2 / P4: a single-phase body beside a two-phase one; the grid column only the single-phase
    particles reach keeps p = 0 (fixed), the run stays finite."""
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 1.0, "-cells", 2, 2, "-basis", "gimp")
    ops.mpm("mode", "implicit")
    ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 2.0, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    ops.mpm("body", 2, "box", 1.0, 0.0, 2.0, 2.0, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    ops.mpm("twoPhase", 1, "-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-perm", 0.01, 0.01)
    ops.mpm("porePressure", 1, "uniform", 5.0)
    ops.mpm("fix", "ymin", "y")
    ops.mpm("gravity", 0.0, -9.81)
    _analysis(ops)
    ops.mpm("analyze", 2, 0.01)
    p = _as_list(ops.mpm("particleData", "p"))
    body = _as_list(ops.mpm("particleData", "body"))
    assert all(math.isfinite(x) for x in p)
    assert all(pp == 0.0 for pp, b in zip(p, body) if b == 2.0)
    stray = 0
    for t in _as_list(ops.mpm("info", "gridNodeTags")):
        if ops.nodeCoord(int(t))[0] > 1.5:
            stray += 1
            assert ops.nodeDisp(int(t), 3) == 0.0
    assert stray > 0


def _restart_deck(ops):
    _deck(ops)
    ops.mpm("porePressure", 1, "uniform", 5.0)
    ops.mpm("gravity", 0.0, -9.81)
    ops.mpm("plate", "ymax", 4.0)
    _analysis(ops)


def test_t9_restart_bit_for_bit(tmp_path):
    """T9: a database round trip mid-consolidation continues bit for bit (pack 5: p, pd, pdd, the
    body's fluid constants, pfix, the plate and its node)."""
    ops = _ops()
    _restart_deck(ops)
    ops.mpm("analyze", 6, 0.01)
    ref = {f: _as_list(ops.mpm("particleData", f)) for f in ("x", "v", "sig", "p", "pd", "pdd")}
    _restart_deck(ops)
    ops.mpm("analyze", 3, 0.01)
    ops.database("File", str(tmp_path / "tpdb"))
    ops.save(3)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.database("File", str(tmp_path / "tpdb"))
    ops.restore(3)
    _analysis(ops)
    ops.mpm("analyze", 3, 0.01)
    for f, r in ref.items():
        assert _as_list(ops.mpm("particleData", f)) == r, f


def test_t12_engines_agree(tmp_path):
    """T12: classic Tcl and Python give the same answer on a small two-phase column (a short T2)."""
    import os
    import subprocess
    from pathlib import Path
    ops = _ops()
    exe = os.environ.get("LADRUNO_TCL_EXE") or str(Path(__file__).resolve().parents[1] / "dist" / "bin" /
                                                   "OpenSees.exe")
    if not os.path.isfile(exe):
        pytest.skip("no OpenSees.exe")
    deck = [
        "model basic -ndm 2 -ndf 2",
        "nDMaterial ElasticIsotropic 1 1.0e4 0.2 2.0",
        "mpm grid 232 -origin 0.0 0.0 -h 1.0 -cells 1 4 -basis gimp",
        "mpm mode implicit",
        "mpm body 1 box 0.0 0.0 1.0 4.0 -ppc 2 -mat 1 -rho 2.0",
        "mpm twoPhase 1 -Kf 4444.0 -poro 0.4 -rhoF 1.0 -perm 1.0e-4 1.0e-4 -stab auto 0.25",
        "mpm fix xmin x", "mpm fix xmax x", "mpm fix ymin y", "mpm pfix ymax 0.0",
        "mpm plate ymax 10.0",
        "constraints Transformation", "numberer RCM", "system UmfPack",
        "test NormDispIncr 1.0e-10 25", "algorithm Newton", "integrator Newmark 0.6 0.3025",
        "analysis Transient",
        "mpm analyze 20 0.02",
    ]
    tcl = tmp_path / "t12.tcl"
    tcl.write_text("\n".join(deck + ["puts \"P [mpm particleData p]\"", "puts \"X [mpm particleData x]\""]) + "\n")
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    proc = subprocess.run([exe, str(tcl)], capture_output=True, text=True, timeout=300, env=env,
                          stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr
    tp = [float(v) for line in out.splitlines() if line.startswith("P ") for v in line.split()[1:]]
    tx = [float(v) for line in out.splitlines() if line.startswith("X ") for v in line.split()[1:]]
    ops.wipe()
    for line in deck:
        t = line.split()
        args = [_num(a) for a in t[1:]]
        getattr(ops, t[0])(*args)
    pp = _as_list(ops.mpm("particleData", "p"))
    px = _as_list(ops.mpm("particleData", "x"))
    assert len(tp) == len(pp) > 0 and len(tx) == len(px), out[-2000:]
    assert max(abs(a - b) for a, b in zip(tp, pp)) <= 1e-12 * max(abs(v) for v in pp)
    assert max(abs(a - b) for a, b in zip(tx, px)) <= 1e-12 * max(abs(v) for v in px)


def _num(a):
    try:
        return int(a)
    except ValueError:
        try:
            return float(a)
        except ValueError:
            return a
