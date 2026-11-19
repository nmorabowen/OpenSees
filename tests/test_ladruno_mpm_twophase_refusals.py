"""WP-232l T8 (part): two-phase refusals and options, mode A (interface P7, §7a Q4).

- `mpm massCutoff` leaves route S's `mpm pressureCutoff` unchanged (ruling c84f71fb3), and
  `mpm info twoPhase` reports it.
- Loud refusals at `mpm analyze`: a symmetric-storage system, mode quasiStatic, HHT, an
  initial-stress offset on a two-phase particle, `mpm antiLocking fbar`, `-dynSeepage on`.
- Route S (mode F) runs the admitted deck and refuses a wall, F-bar and fix 2's term.
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


def _refused(ops, capfd, needle, call=None):
    """The call raises AND the refusal names its rule: `needle` is in the stderr text (review of
    record: no bare raises)."""
    capfd.readouterr()
    with pytest.raises(Exception):
        if call is None:
            ops.mpm("analyze", 1, 0.01)
        else:
            call()
    err = capfd.readouterr().err
    assert needle in err, err[-800:]


REFUSAL_TEXT = {"symmetric": "symmetric-storage", "static": "mode quasiStatic is refused",
                "hht": "Newmark only", "sig0": "initial-stress offset", "fbar": "antiLocking fbar",
                "masscut": "massCutoff"}


def test_pressure_cutoff_is_its_own_knob():
    ops = _ops()
    _deck(ops)
    assert _as_list(ops.mpm("info", "twoPhase", 1))[11] == 5e-2
    ops.mpm("massCutoff", 1e-10)
    assert _as_list(ops.mpm("info", "twoPhase", 1))[11] == 5e-2
    ops.mpm("pressureCutoff", 0.1)
    assert _as_list(ops.mpm("info", "twoPhase", 1))[11] == 0.1
    ops.mpm("massCutoff", 0.3)
    assert _as_list(ops.mpm("info", "twoPhase", 1))[11] == 0.1


def test_admitted_deck_runs():
    ops = _ops()
    _deck(ops)
    _analysis(ops)
    ops.mpm("analyze", 1, 0.01)
    p = _as_list(ops.mpm("particleData", "p"))
    assert len(p) == 16 and all(x == x for x in p)


@pytest.mark.parametrize("why", ["symmetric", "static", "hht", "sig0", "fbar", "masscut"])
def test_refusals(why, capfd):
    """... "masscut": Fable audit finding 4, mode A refuses massCutoff > pressureCutoff."""
    ops = _ops()
    extra = {"fbar": [("antiLocking", "fbar")], "masscut": [("massCutoff", 0.2)]}.get(why, [])
    _deck(ops, mode="quasiStatic" if why == "static" else "implicit", stress=(why == "sig0"),
          extra=extra)
    if why == "symmetric":
        _analysis(ops, system="ProfileSPD")
    elif why == "static":
        _analysis(ops, static=True)
    elif why == "hht":
        _analysis(ops, integrator=("HHT", 0.9))
    else:
        _analysis(ops)
    _refused(ops, capfd, REFUSAL_TEXT[why])


def test_dynseepage_on_refused(capfd):
    ops = _ops()
    _deck(ops)
    _refused(ops, capfd, "dynSeepage on is refused",
             lambda: ops.mpm("twoPhase", 1, "-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-perm", 0.01,
                             0.01, "-dynSeepage", "on"))


def test_route_s_admitted_deck_runs():
    """Route S (mode F) runs the admitted two-phase deck (gimp, pfix top, base fixed)."""
    ops = _ops()
    _deck(ops, mode="explicit")
    ops.mpm("analyze", 2, 1.0e-3)
    p = _as_list(ops.mpm("particleData", "p"))
    assert len(p) == 16 and all(x == x for x in p)
    assert _as_list(ops.mpm("info", "twoPhase", 1))[0] == 2.0            # route S


@pytest.mark.parametrize("why", ["wall", "fbar", "fix2", "sig0"])
def test_route_s_refusals(why, capfd):
    """Route S refuses what it does not carry (P7): a wall, F-bar, and fix 2's affine term (a
    B-spline impervious face with a normal fluid drive)."""
    ops = _ops()
    if why == "fix2":
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
        ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 0.5, "-cells", 2, 4, "-basis", "bspline2")
        ops.mpm("mode", "explicit")
        ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 2.0, "-ppc", 2, "-mat", 1, "-rho", 2.0)
        ops.mpm("twoPhase", 1, "-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-perm", 0.01, 0.01)
        ops.mpm("gravity", 0.0, -9.81)
        ops.mpm("fix", "ymin", "y")
        ops.mpm("pfix", "ymax", 0.0)
    elif why == "sig0":                                  # mode F: a sigma0 offset on a two-phase particle
        _deck(ops, mode="explicit", stress=True)
    else:
        _deck(ops, mode="explicit",
              extra=[("wall", "xmax", 0.0)] if why == "wall" else [("antiLocking", "fbar")])
    needle = {"wall": "`mpm wall` is refused", "fbar": "antiLocking fbar",
              "fix2": "fix 2's affine image term", "sig0": "initial-stress offset"}[why]
    _refused(ops, capfd, needle, lambda: ops.mpm("analyze", 1, 1.0e-3))


def test_generalized_alpha_refused(capfd):
    ops = _ops()
    _deck(ops)
    _analysis(ops, integrator=("GeneralizedAlpha", 0.9, 0.8))
    _refused(ops, capfd, "Newmark only")


@pytest.mark.parametrize("args,needle", [
    (["-poro", 0.4, "-rhoF", 1.0, "-perm", 0.01, 0.01], "missing required argument"),          # no -Kf
    (["-Kf", 4444.0, "-poro", 0.6, "-rhoF", 1.0, "-perm", 0.01, 0.01, "-alpha", 0.5], "alpha"),   # n > alpha
    (["-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-perm", 0.01, 0.01, "-bogus", 1.0], "unknown flag"),
    (["-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-permH", 0.01, 0.01], "-permH and -gammaW come together"),
])
def test_twophase_parser_refusals(args, needle, capfd):
    ops = _ops()
    _deck(ops)
    _refused(ops, capfd, needle, lambda: ops.mpm("twoPhase", 2, *args))


def test_porous_material_refused_at_particles(capfd):
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    ops.nDMaterial("FluidSolidPorous", 2, 2, 1, 2.2e6, 101.0)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 1.0, "-cells", 2, 2)
    _refused(ops, capfd, "Porous materials (PDMY, PIMY, FluidSolidPorous) are not admitted",
             lambda: ops.mpm("body", 1, "box", 0.0, 0.0, 2.0, 2.0, "-ppc", 2, "-mat", 2, "-rho", 2.0))


def test_bspline_pfix_nonzero_and_traction_refused(capfd):
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 1.0, "-cells", 2, 2, "-basis", "bspline2")
    ops.mpm("pfix", "ymax", 0.0)                     # drained: admitted
    _refused(ops, capfd, "p0 != 0 on a B-spline face is refused", lambda: ops.mpm("pfix", "xmax", 3.0))
    _refused(ops, capfd, "is GIMP / linear only in phase 6",
             lambda: ops.mpm("traction", "t", "ymax", 0.0, 2.0, 1.0, 2.0, 0.0, -1.0))


def test_pack_below_5_refused_with_two_phase(capfd):
    ops = _ops()
    _deck(ops)
    _refused(ops, capfd, "two-phase", lambda: ops.mpm("debug", "packVersion", 4))


def test_info_reports_route_constants_and_stab():
    ops = _ops()
    _deck(ops)
    info = _as_list(ops.mpm("info", "twoPhase", 1))
    # [route, Kf, n, rhoF, k_x, k_y, alpha, Ks, stabMode, a0, alpha_st(h), pressureCutoff]
    assert info[:10] == [1.0, 4444.0, 0.4, 1.0, 0.01, 0.01, 1.0, -1.0, 1.0, 0.25]
    G = 1.0e4 / 2.4
    lam = 1.0e4 * 0.2 / (1.2 * 0.6)
    assert info[10] == pytest.approx(0.25 / (lam + 2 * G / 3 + 4 * G / 3), rel=1e-12)


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


def test_rank_deficient_projection_refused(capfd):
    """P5a's singularity guard (b): bspline2 with one particle per cell -- fewer particles than
    p-active nodes, so the consistent projection matrix is singular; the step is refused loudly
    and nothing moves."""
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 1.0, "-cells", 4, 4, "-basis", "bspline2")
    ops.mpm("mode", "implicit")
    ops.mpm("body", 1, "box", 0.0, 0.0, 4.0, 4.0, "-ppc", 1, "-mat", 1, "-rho", 2.0)
    ops.mpm("twoPhase", 1, "-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-perm", 0.01, 0.01)
    ops.mpm("porePressure", 1, "uniform", 5.0)
    ops.mpm("fix", "ymin", "y")
    _analysis(ops)
    x0 = _as_list(ops.mpm("particleData", "x"))
    _refused(ops, capfd, "consistent pressure projection")
    assert _as_list(ops.mpm("particleData", "x")) == x0
    assert all(v == 5.0 for v in _as_list(ops.mpm("particleData", "p")))


def test_sealed_column_keeps_p_bit_for_bit():
    """P5's invariant, a guard on P5a: a motionless, sealed, kbar = 0 column keeps p, pd, pdd bit
    for bit (no gravity, no load, no drained face, a non-uniform p)."""
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 1.0, "-cells", 2, 4, "-basis", "gimp")
    ops.mpm("mode", "implicit")
    ops.mpm("body", 1, "box", 0.0, 0.0, 2.0, 4.0, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    ops.mpm("twoPhase", 1, "-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-perm", 0.0, 0.0)
    n = len(_as_list(ops.mpm("particleData", "m")))
    p0 = [5.0 + 0.37 * math.sin(1.3 * k) for k in range(n)]
    ops.mpm("particlePressure", *(p0 + [0.0] * n + [0.0] * n))
    # motionless: every grid node held at zero velocity (a homogeneous gridvel over the box)
    ops.mpm("gridvel", "all", "box", -1.5, -1.5, 3.5, 5.5, "v", 0.0, 0.0)
    _analysis(ops)
    ops.mpm("analyze", 5, 0.01)
    assert _as_list(ops.mpm("particleData", "p")) == p0
    assert all(v == 0.0 for v in _as_list(ops.mpm("particleData", "pd")))


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


def test_drained_face_p_weights_reproduce_linear_after_u_drop(tmp_path):
    """Ruling bb1421559 (unit level): the row beyond a drained (pfix) B-spline face is dropped on
    the u side (mass below `mpm massCutoff`); the u weights renormalise, but the p weights keep
    fix 3's rule on the raw weights (the folded-away ghost counts as retained), so at the top
    particle they still reproduce p = (H - y), the linear field odd about the face, to 1e-14."""
    import json
    ops = _ops()
    H, h = 10.0, 0.5
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", h, "-cells", 2, 20, "-basis", "bspline2")
    ops.mpm("mode", "implicit")
    ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, H, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    ops.mpm("twoPhase", 1, "-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-perm", 1.0, 1.0, "-stab", "off")
    for face, c in (("xmin", "x"), ("xmax", "x"), ("ymin", "y")):
        ops.mpm("fix", face, c)
    ops.mpm("pfix", "ymax", 0.0)
    ops.mpm("massCutoff", 0.2)              # drops the y = H + h row (its mass ~0.05 m_max)
    ops.mpm("pressureCutoff", 0.25)         # finding 4: massCutoff <= pressureCutoff in mode A
    _analysis(ops)
    xy = _as_list(ops.mpm("particleData", "x"))
    top = max(range(len(xy) // 2), key=lambda k: (xy[2 * k + 1], -xy[2 * k]))
    y = xy[2 * top + 1]
    xi = (H + h - y) / h
    assert 0.0 < 0.5 * (1.5 - xi) ** 2       # the raw weight of the dropped row is not zero
    path = str(tmp_path / "top.json")
    ops.mpm("debug", "particleBlocks", top, path, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    d = json.loads(open(path).read())
    Y = [c[1] for c in d["coords"]]
    assert max(Y) <= H + 1e-9, "the row beyond the drained face was not u-dropped"
    assert abs(sum(d["S"]) - 1.0) <= 1e-14  # the u weights renormalised
    lin = [H - v for v in Y]
    f = sum(s * v for s, v in zip(d["Sp"], lin))
    gy = sum(g[1] * v for g, v in zip(d["dSp"], lin))
    assert abs(f - (H - y)) <= 1e-14
    assert abs(gy + 1.0) <= 1e-13


def test_plate_pfix_face_weights_sum_to_one(tmp_path):
    """Fable audit finding 1 (wp-briefs "232l"): beside a face with a u fold (`plate`) and `pfix`,
    the ghost's W folds onto the mirror only -- sum W over the lattice entries is 1 to 1e-14 at
    every particle within reach of the face (`emptyNodes exclude`, `massCutoff 0`: no
    renormalisation hides it). The dangling else gave 1.03125."""
    import json
    ops = _ops()
    H, h = 10.0, 1.0
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", h, "-cells", 2, 10, "-basis", "bspline2")
    ops.mpm("mode", "implicit")
    ops.mpm("body", 1, "box", 0.0, 0.0, 2.0, H, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    ops.mpm("twoPhase", 1, "-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-perm", 1e-4, 1e-4, "-stab", "off")
    for face, c in (("xmin", "x"), ("xmax", "x"), ("ymin", "y")):
        ops.mpm("fix", face, c)
    ops.mpm("pfix", "ymax", 0.0)
    ops.mpm("plate", "ymax", 20.0)
    ops.mpm("emptyNodes", "exclude")
    ops.mpm("massCutoff", 0.0)
    _analysis(ops)
    xy = _as_list(ops.mpm("particleData", "x"))
    checked = 0
    for p in range(len(xy) // 2):
        if xy[2 * p + 1] < H - 1.5 * h:
            continue
        path = str(tmp_path / ("p%d.json" % p))
        ops.mpm("debug", "particleBlocks", p, path, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        d = json.loads(open(path).read())
        grid = set(int(t) for t in _as_list(ops.mpm("info", "gridNodeTags")))   # after the dump's build
        sw = sum(s for t, s in zip(d["nodes"], d["S"]) if int(t) in grid)
        assert abs(sw - 1.0) <= 1e-14, (p, xy[2 * p + 1], sw)
        checked += 1
    assert checked > 0



def _column(ops, ramp=False, newton_fail=False):
    """A small bspline2 hydrostatic column (the fix-2 patch deck, H = 4): pi != 0 under gravity."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 0.5, "-cells", 2, 8, "-basis", "bspline2")
    ops.mpm("mode", "implicit")
    ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 4.0, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    ops.mpm("twoPhase", 1, "-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-perm", 1.0, 1.0, "-stab", "off")
    ops.mpm("gravity", 0.0, -9.81)
    for face, c in (("xmin", "x"), ("xmax", "x"), ("ymin", "y")):
        ops.mpm("fix", face, c)
    ops.mpm("pfix", "ymax", 0.0)
    if ramp:
        ops.timeSeries("Path", 77, "-time", 0.0, 5.0, 100.0, "-values", 0.0, 1.0, 1.0)
        ops.mpm("timeSeries", 77)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    if newton_fail:
        ops.test("NormDispIncr", 1e-300, 1)
        ops.mpm("cut", 2.0, 0)
    else:
        ops.test("NormDispIncr", 1e-10, 25)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")


def _pi_perp_max(ops):
    return _as_list(ops.mpm("info", "twoPhase", 1))[-1]


def test_a_prime_commits_on_acceptance_only(tmp_path):
    """Fable audit finding 2: A' (pi_perp) and pi_c are committed only when a step is accepted.
    Block dumps change nothing (p bit for bit, max|pi_perp| 0); dump -> analyze gives the same p
    and pi_perp as a fresh analyze; reset -> analyze gives them again (one pi_perp)."""
    ops = _ops()
    _column(ops)
    ops.mpm("analyze", 1, 1.0)
    p_ref, pp_ref = _as_list(ops.mpm("particleData", "p")), _pi_perp_max(ops)
    assert pp_ref > 0.0
    _column(ops)
    p0 = _as_list(ops.mpm("particleData", "p"))
    for k in (0, 5):
        ops.mpm("debug", "particleBlocks", k, str(tmp_path / ("d%d.json" % k)), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    assert _as_list(ops.mpm("particleData", "p")) == p0
    assert _pi_perp_max(ops) == 0.0
    ops.mpm("analyze", 1, 1.0)
    p1 = _as_list(ops.mpm("particleData", "p"))
    assert max(abs(a - b) for a, b in zip(p1, p_ref)) <= 1e-12 * 9.81 * 4.0
    assert abs(_pi_perp_max(ops) - pp_ref) <= 1e-12 * 9.81 * 4.0
    ops.mpm("reset")
    ops.mpm("analyze", 1, 1.0)
    p2 = _as_list(ops.mpm("particleData", "p"))
    assert max(abs(a - b) for a, b in zip(p2, p_ref)) <= 1e-12 * 9.81 * 4.0
    assert abs(_pi_perp_max(ops) - pp_ref) <= 1e-12 * 9.81 * 4.0


def test_refused_first_step_leaves_a_prime_unset():
    """Fable audit finding 2: a first step that fails (Newton, no cut allowed) leaves p and the A'
    flag untouched; the next accepted step applies pi_perp once, as a fresh run does."""
    ops = _ops()
    _column(ops)
    ops.mpm("analyze", 1, 1.0)
    p_ref, pp_ref = _as_list(ops.mpm("particleData", "p")), _pi_perp_max(ops)
    _column(ops, newton_fail=True)
    p0 = _as_list(ops.mpm("particleData", "p"))
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 1.0)
    assert _as_list(ops.mpm("particleData", "p")) == p0
    assert _pi_perp_max(ops) == 0.0
    ops.test("NormDispIncr", 1e-10, 25)
    ops.mpm("analyze", 1, 1.0)
    p1 = _as_list(ops.mpm("particleData", "p"))
    assert max(abs(a - b) for a, b in zip(p1, p_ref)) <= 1e-12 * 9.81 * 4.0
    assert abs(_pi_perp_max(ops) - pp_ref) <= 1e-12 * 9.81 * 4.0


def test_ramped_drive_needs_no_a_prime():
    """Fable audit finding 3: pi is evaluated at the drive of the step being taken; A' compensates
    pi(0) at the epoch's initial drive, which is 0 for a ramp from 0: max|pi_perp| = 0."""
    ops = _ops()
    _column(ops, ramp=True)
    ops.mpm("analyze", 2, 1.0)
    assert _pi_perp_max(ops) == 0.0
