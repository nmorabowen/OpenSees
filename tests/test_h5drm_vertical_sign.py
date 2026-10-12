# ============================================================================
# test_h5drm_vertical_sign.py  (WP-186)
# ----------------------------------------------------------------------------
# End-to-end sign test of the H5DRM load pattern's vertical component.
#
# Convention under test (H5DRMLoadPattern.cpp, WP-186): an .h5drm file is in
# ShakerMaker's frame -- station xyz = (N, E, depth) in km, data rows
# (E, N, vertical) with the vertical positive DOWN. T and x0 act on the station
# coordinates only; the rows enter the model as u_X = +d0, u_Y = +d1,
# u_Z = -d2. With the decks' T = (0 1 0)(1 0 0)(0 0 -1) the model is z-up, so
# a z-down file moves the model the right way, and a file written "vertical up"
# (the Paper 1 Legacy convention, built for the no-flip binary) moves it the
# wrong way.
#
# Exact fixture: the free field is a rigid translation u0(t) = r(t) * U with
# r = t / tend, zero acceleration, solved statically. A rigid translation is an
# exact static free field (zero strain, traction-free surface), so the DRM
# effective forces must reproduce it to round-off in the interior (b nodes and
# inside) and leave zero scattered field outside (e nodes and buffer). Units and
# the pattern line mirror the Paper 1 decks: model in mm, file in m / km,
#   pattern H5DRM tag file 1000 1e6 2500 1  0 1 0 1 0 0 0 0 -1  22000 15500 0
#
# Physical motion: surface UP 1 mm, East 2 mm, North 3 mm (all scaled by r).
#   z-down file (row 2 = -1 mm)  -> surface u_Z = +1 mm * r   (correct)
#   up file     (row 2 = +1 mm)  -> surface u_Z = -1 mm * r   (inverted)
# Horizontal u_X = +2 mm * r (East), u_Y = +3 mm * r (North) in both.
# At t = tend the pattern holds the final displacement sample: the hold branch must
# carry the same flip (no sign jump at tend).
#
# Runs through openseespy (worktree dist/bin/opensees.pyd) and, when present,
# the classic Tcl OpenSees.exe on the same file (the parser the decks use).
# ============================================================================
import os
import subprocess
import sys

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_DIST = os.path.normpath(os.path.join(_HERE, "..", "dist", "bin"))

h5py = pytest.importorskip("h5py")
if not os.path.exists(os.path.join(_DIST, "opensees.pyd")):
    pytest.skip("worktree dist/bin/opensees.pyd not built", allow_module_level=True)
os.add_dll_directory(_DIST)
sys.path.insert(0, _DIST)
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
import opensees as ops  # noqa: E402

# --- fixture geometry (ShakerMaker frame, km) --------------------------------
H_KM = 0.01                      # 10 m station spacing
NI = 2                           # stations at -2h..2h in N and E
ND = 3                           # depths 0..3h (z positive down)
DT = 0.01
NT = 101                         # samples; tend = 1.0 s
TEND = 1.0
# physical free field at full ramp, metres: East, North, UP
U_E, U_N, U_UP = 2.0e-3, 3.0e-3, 1.0e-3

# --- deck mapping (as the rev2_m36 Tier 4 decks) ------------------------------
FACTOR = 1000.0                  # m -> mm (displacements)
CRD = 1.0e6                      # km -> mm (coordinates)
DTOL = 2500.0                    # mm
T = (0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0)
X0 = (22000.0, 15500.0, 0.0)
HMM = H_KM * CRD                 # 10 000 mm


def _stations():
    out = []
    for i in range(-NI, NI + 1):          # x = N
        for j in range(-NI, NI + 1):      # y = E
            for k in range(0, ND + 1):    # z = depth
                out.append((i, j, k))
    return out


def _is_internal(i, j, k):
    return abs(i) <= NI - 1 and abs(j) <= NI - 1 and k <= ND - 1


def write_h5drm(path, vertical):
    """vertical='down': ShakerMaker native (row 2 positive down);
    vertical='up': Legacy as-run convention (row 2 positive up)."""
    st = _stations()
    n = len(st)
    xyz = np.array([(i * H_KM, j * H_KM, k * H_KM) for i, j, k in st], dtype=float)
    internal = np.array([_is_internal(*s) for s in st], dtype=bool)
    r = np.arange(NT) * DT / TEND
    zrow = -U_UP if vertical == "down" else +U_UP
    dis = np.zeros((3 * n, NT))
    for p in range(n):
        dis[3 * p + 0] = U_E * r
        dis[3 * p + 1] = U_N * r
        dis[3 * p + 2] = zrow * r
    zeros = np.zeros_like(dis)
    with h5py.File(path, "w") as f:
        d = f.create_group("DRM_Data")
        d["xyz"] = xyz
        d["internal"] = internal
        d["data_location"] = (np.arange(n) * 3).astype(np.int32)
        d["displacement"] = dis
        d["velocity"] = zeros
        d["acceleration"] = zeros
        m = f.create_group("DRM_Metadata")
        m["dt"] = DT
        m["tstart"] = 0.0
        m["tend"] = TEND
        m["drmbox_x0"] = np.array([0.0, 0.0, 0.0])
        for kname, v in (("drmbox_xmin", -NI * H_KM), ("drmbox_xmax", NI * H_KM),
                         ("drmbox_ymin", -NI * H_KM), ("drmbox_ymax", NI * H_KM),
                         ("drmbox_zmin", 0.0), ("drmbox_zmax", ND * H_KM)):
            m[kname] = float(v)
        m["h"] = np.array([H_KM, H_KM, H_KM])


def model_grid():
    """Model nodes (mm, z up): the station grid plus one buffer layer (|i| or
    |j| = 3, depth 4h) that is fixed and carries no station. Returns
    {(i, j, k): (tag, X, Y, Z)} with i = N index, j = E index, k = depth index."""
    nodes = {}
    tag = 1
    for i in range(-NI - 1, NI + 2):
        for j in range(-NI - 1, NI + 2):
            for k in range(0, ND + 2):
                X = j * HMM + X0[0]          # East
                Y = i * HMM + X0[1]          # North
                Z = -k * HMM                 # up
                nodes[(i, j, k)] = (tag, X, Y, Z)
                tag += 1
    return nodes


def bricks(nodes):
    """stdBrick connectivity over every grid cell (bottom face CCW, then top)."""
    out = []
    for i in range(-NI - 1, NI + 1):
        for j in range(-NI - 1, NI + 1):
            for k in range(0, ND + 1):
                # model X ~ j, Y ~ i, Z ~ -k ; bottom = deeper (k+1)
                def t(di, dj, dk):
                    return nodes[(i + di, j + dj, k + dk)][0]
                out.append([t(0, 0, 1), t(0, 1, 1), t(1, 1, 1), t(1, 0, 1),
                            t(0, 0, 0), t(0, 1, 0), t(1, 1, 0), t(1, 0, 0)])
    return out


def is_buffer(i, j, k):
    return abs(i) == NI + 1 or abs(j) == NI + 1 or k == ND + 1


def run_openseespy(h5):
    nodes = model_grid()
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for (i, j, k), (tag, X, Y, Z) in nodes.items():
        ops.node(tag, X, Y, Z)
        if is_buffer(i, j, k):
            ops.fix(tag, 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, 200.0, 0.3, 2.0e-9)
    for e, ns in enumerate(bricks(nodes), start=1):
        ops.element("stdBrick", e, *ns, 1)
    ops.pattern("H5DRM", 1, h5, FACTOR, CRD, DTOL, 1, *T, *X0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-10, 5)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 0.5)                 # static: no inertia, exact rigid field
    ops.analysis("Static")
    out = {}
    assert ops.analyze(1) == 0                           # t = 0.5, inside the record
    out["mid"] = {key: ops.nodeDisp(nodes[key][0]) for key in nodes}
    assert ops.analyze(1) == 0                           # t = 1.0 = tend exactly: the hold branch
    out["hold"] = {key: ops.nodeDisp(nodes[key][0]) for key in nodes}
    out["t"] = ops.getTime()
    ops.wipe()
    return nodes, out


def expected(vertical, r):
    zs = +1.0 if vertical == "down" else -1.0
    return np.array([U_E, U_N, zs * U_UP]) * FACTOR * r


@pytest.mark.parametrize("vertical", ["down", "up"])
def test_h5drm_vertical_sign_openseespy(tmp_path, vertical):
    h5 = str(tmp_path / f"rigid_{vertical}.h5drm")
    write_h5drm(h5, vertical)
    nodes, out = run_openseespy(h5)
    assert out["t"] == TEND
    surf = (0, 0, 0)                                     # interior surface centre
    u_mid = np.array(out["mid"][surf])
    exp_mid = expected(vertical, 0.5)
    np.testing.assert_allclose(u_mid, exp_mid, rtol=1e-6, atol=1e-9)
    if vertical == "down":
        assert u_mid[2] > 0.0, "z-down file must move the z-up model UP"
    else:
        assert u_mid[2] < 0.0, "an up-convention file is inverted by the flip"
    # every interior node (b layer included) moves rigidly with the free field
    for key, u in out["mid"].items():
        if _is_internal(*key) and not is_buffer(*key):
            np.testing.assert_allclose(u, exp_mid, rtol=1e-6, atol=1e-9)
    # the exterior (e layer + buffer) carries no scattered field
    for key, u in out["mid"].items():
        if not _is_internal(*key):
            assert np.max(np.abs(u)) < 1e-6 * np.max(np.abs(exp_mid)), key
    # t = tend: drm_direct_read holds sample i_final = N_timesteps - 1 (r = 0.99) with the same
    # flip. (Past tend, CalculateBoundaryForces returns zero DRM forces: t > tend is not reached.)
    r_hold = (round(TEND / DT) - 1) * DT / TEND
    u_hold = np.array(out["hold"][surf])
    np.testing.assert_allclose(u_hold, expected(vertical, r_hold), rtol=1e-6, atol=1e-9)


_EXE = os.path.join(_DIST, "OpenSees.exe")


@pytest.mark.skipif(not os.path.exists(_EXE), reason="dist/bin/OpenSees.exe not built")
@pytest.mark.parametrize("vertical", ["down", "up"])
def test_h5drm_vertical_sign_tcl(tmp_path, vertical):
    """Same model through classic OpenSees.exe and the deck's pattern syntax."""
    h5 = (tmp_path / f"rigid_{vertical}.h5drm")
    write_h5drm(str(h5), vertical)
    nodes = model_grid()
    lines = ["wipe", "model basic -ndm 3 -ndf 3"]
    for (i, j, k), (tag, X, Y, Z) in nodes.items():
        lines.append(f"node {tag} {X!r} {Y!r} {Z!r}")
        if is_buffer(i, j, k):
            lines.append(f"fix {tag} 1 1 1")
    lines.append("nDMaterial ElasticIsotropic 1 200.0 0.3 2.0e-9")
    for e, ns in enumerate(bricks(nodes), start=1):
        lines.append(f"element stdBrick {e} {' '.join(map(str, ns))} 1")
    tt = " ".join(repr(v) for v in T)
    xx = " ".join(repr(v) for v in X0)
    lines.append(f'pattern H5DRM 1 "{h5.as_posix()}" {FACTOR} {CRD} {DTOL} 1  {tt}  {xx}')
    s = nodes[(0, 0, 0)][0]
    out = (tmp_path / "u.txt").as_posix()
    lines += ["constraints Transformation", "numberer RCM", "system UmfPack",
              "test NormDispIncr 1e-10 5", "algorithm Linear",
              "integrator LoadControl 0.5", "analysis Static",
              "if {[analyze 1] != 0} {exit 2}",
              f'set f [open "{out}" w]; puts $f [nodeDisp {s}]; close $f', "exit 0"]
    deck = tmp_path / "deck.tcl"
    deck.write_text("\n".join(lines) + "\n")
    p = subprocess.run([_EXE, str(deck)], cwd=str(tmp_path), capture_output=True, text=True,
                       timeout=300)
    assert p.returncode == 0, p.stdout[-2000:] + p.stderr[-2000:]
    u = np.array([float(v) for v in open(out).read().split()])
    np.testing.assert_allclose(u, expected(vertical, 0.5), rtol=1e-6, atol=1e-9)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
