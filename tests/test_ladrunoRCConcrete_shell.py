"""LadrunoRCConcrete in its REAL target element — ASDShellQ4 + LayeredShellFiberSection
(classTag 33015, PlateFiber view via getCopy) — Zone-A shell-integration battery.

Phases 1/2a/2b.1 were all verified on a 3D stdBrick (the full 6-component material). This
file closes the gap the adversarial reviews flagged: the material is a SHELL material, and
the PlateFiber view (order-5, with the guarded sigma_33=0 inner Newton) had never been
exercised inside an actual ASDShellQ4 + LayeredShell section end-to-end.

Driver: a single flat ASDShellQ4 unit quad whose every membrane DOF is prescribed (Penalty)
to realize a HOMOGENEOUS membrane strain; we read the section stress resultants
(eleResponse 'stresses' -> [Nxx,Nyy,Nxy,Mxx,Myy,Mxy,Vxz,Vyz] per Gauss point). With a uniform
material across layers and a uniform membrane strain every layer is identical, so the resultant
is N = (layer stress) * h.

Validated here, end-to-end through the shell:
  * membrane tension cracks + softens (Nxx peaks ~ ft*h at eps_xx ~ ft/E, then drops);
  * cyclic membrane shear engages the fixed-crack aggregate interlock: Nxy saturates at the
    MCFT bound v_ci,max * h on BOTH signs (PlateFiber x sigma_33-condensation x cyclic);
  * the interlock is load-bearing in the shell (ON caps Nxy below the OFF value).

Plan: Ladruno_implementation/19_ladruno_rc_shell_adr.md (Phase 2b.2).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# same curved backbone as the brick battery
E, NU, KC = 30000.0, 0.2, 2.0 / 3.0
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0,   30.0,   5.0]
CD = [0.0, 0.0,    0.25,   1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010]
TS = [0.0, 3.0,    0.5]
TD = [0.0, 0.0,    1.0 - 0.5 / 5.0]
FT = max(TS)                       # tensile peak
H = 0.1                            # shell thickness
NLAYERS = 4
_QUAD = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}


def _mat(tag, interlock=False, cyclic=False, crack_spacing=0.0, agg=16.0):
    args = ["LadrunoRCConcrete", tag, E, NU,
            "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
            "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC]
    if interlock:
        args += ["-interlock", "-agg", agg]
        if crack_spacing > 0.0:
            args += ["-crackSpacing", crack_spacing]
        if cyclic:
            args += ["-cyclic"]
    ops.nDMaterial(*args)


def _build_shell(mat_fn):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y) in _QUAD.items():
        ops.node(t, x, y, 0.0)
    mat_fn(1)
    layers = []
    for _ in range(NLAYERS):
        layers += [1, H / NLAYERS]
    ops.section("LayeredShell", 1, NLAYERS, *layers)
    ops.element("ASDShellQ4", 1, 1, 2, 3, 4, 1)


def _membrane_resultants():
    ops.eleResponse(1, "forces")                          # prime lazy state
    st = list(ops.eleResponse(1, "stresses"))             # 4 GP x 8 resultants
    return st[0], st[1], st[2]                             # Nxx, Nyy, Nxy at GP0


def _prescribe_membrane(exx, gxy):
    """u_x = exx*X + gxy*Y? No -- use the non-symmetric simple shear u_y=gxy*X (disjoint from
    the exx*X tension on dof-1) so a per-stage Path can ramp them independently; fix u_z + the
    three rotation DOFs. (engineering gxy = du_y/dx = gxy.)"""
    for t, (x, y) in _QUAD.items():
        ops.sp(t, 1, exx * x)     # u_x = exx*X
        ops.sp(t, 2, gxy * x)     # u_y = gxy*X  -> engineering shear
        ops.sp(t, 3, 0.0)         # u_z
        ops.sp(t, 4, 0.0); ops.sp(t, 5, 0.0); ops.sp(t, 6, 0.0)   # rotations


def _vci_max(sqrt_fc, w, a_g):
    return 0.18 * sqrt_fc / (0.31 + 24.0 * w / (a_g + 16.0))


# --------------------------------------------------------------------------
@pytest.mark.t1
def test_shell_membrane_tension_cracks_and_softens():
    """Homogeneous membrane tension on the layered RC shell: the membrane resultant Nxx must
    peak near ft*h at the cracking strain (~ft/E) and then SOFTEN -- the full
    ASDShellQ4 -> LayeredShell -> LadrunoRCConcrete(PlateFiber) path under load."""
    _build_shell(lambda t: _mat(t, interlock=False))
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    exx_max, nsteps = 6.0e-4, 60
    for t, (x, y) in _QUAD.items():
        ops.sp(t, 1, exx_max * x)
        ops.sp(t, 2, 0.0); ops.sp(t, 3, 0.0)
        ops.sp(t, 4, 0.0); ops.sp(t, 5, 0.0); ops.sp(t, 6, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e12, 1.0e12)
    ops.test("NormDispIncr", 1.0e-7, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    nxx = []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0, "shell tension analyze failed"
        nxx.append(_membrane_resultants()[0])
    peak = max(nxx)
    assert peak == pytest.approx(FT * H, rel=0.10), f"membrane peak {peak} != ft*h {FT * H}"
    assert nxx[-1] < 0.7 * peak, f"no tension softening: end {nxx[-1]} vs peak {peak}"
    assert nxx[-1] > 0.0, "membrane tension went non-physical (negative)"


def _shell_tension_then_cyclic_shear(mat_fn, exx, gamma, n1=40, nleg=70, ts2_vals=(0, 0, 1, -1, 1)):
    """Two-+-stage membrane path on the shell: stage 1 tension to crack, then shear legs.
    Tension on dof-1 (exx*X), shear on dof-2 (gamma*X). Returns [(gxy, Nxy), ...] over the
    shear legs. Path padded with a hold node (see LEDGER_quirks: Path returns 0 out-of-range)."""
    nst = len(ts2_vals) - 1
    times = [float(i) for i in range(nst + 2)]            # +1 pad node
    ts1 = [0.0] + [1.0] * (nst) + [1.0]                   # tension: ramp in stage 0, hold
    ts2 = list(ts2_vals) + [ts2_vals[-1]]
    _build_shell(mat_fn)
    ops.timeSeries("Path", 1, "-time", *times, "-values", *ts1)
    ops.timeSeries("Path", 2, "-time", *times, "-values", *ts2)
    ops.pattern("Plain", 1, 1)
    for t, (x, y) in _QUAD.items():
        ops.sp(t, 1, exx * x)
        ops.sp(t, 3, 0.0); ops.sp(t, 4, 0.0); ops.sp(t, 5, 0.0); ops.sp(t, 6, 0.0)
    ops.pattern("Plain", 2, 2)
    for t, (x, y) in _QUAD.items():
        ops.sp(t, 2, gamma * x)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e12, 1.0e12)
    ops.test("NormDispIncr", 1.0e-7, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / n1)
    ops.analysis("Static")
    for _ in range(n1):                                   # stage 0: tension to crack
        assert ops.analyze(1) == 0, "shell tension stage failed"
    out = []
    ops.integrator("LoadControl", 1.0 / nleg)
    for stage in range(1, nst):                           # shear legs
        for _ in range(nleg):
            assert ops.analyze(1) == 0, f"shell shear stage {stage} failed"
            out.append((ops.nodeDisp(2, 2), _membrane_resultants()[2]))
    return out


@pytest.mark.t1
def test_shell_cyclic_membrane_shear_caps_both_signs():
    """The headline shell gate: tension cracks the layers, then cyclic membrane shear engages
    the fixed-crack aggregate interlock IN THE ASDShellQ4 element -- the shear resultant Nxy
    saturates at the MCFT bound +/- v_ci,max * h on both half-cycles (PlateFiber view x the
    sigma_33 condensation Newton x the cyclic friction-slip law, end-to-end)."""
    import numpy as np
    exx, gamma, sth, a_g = 5.0e-4, 2.0e-3, 50.0, 16.0
    pts = _shell_tension_then_cyclic_shear(
        lambda t: _mat(t, interlock=True, cyclic=True, crack_spacing=sth, agg=a_g),
        exx, gamma, ts2_vals=(0, 0, 1, -1, 1))
    vci_h = _vci_max(np.sqrt(max(CS)), exx * sth, a_g) * H
    nxy = [n for _, n in pts]
    assert max(nxy) == pytest.approx(vci_h, rel=1.5e-2), f"+Nxy peak {max(nxy)} != v_ci,max*h {vci_h}"
    assert min(nxy) == pytest.approx(-vci_h, rel=1.5e-2), f"-Nxy peak {min(nxy)} != -v_ci,max*h {vci_h}"


@pytest.mark.t1
def test_shell_interlock_is_load_bearing():
    """Ablation in the shell: at large post-crack membrane shear, |Nxy| WITH interlock (capped
    at v_ci,max*h) must be strictly below the interlock-off membrane shear."""
    exx, gamma, sth = 5.0e-4, 3.0e-3, 50.0
    on = _shell_tension_then_cyclic_shear(
        lambda t: _mat(t, interlock=True, cyclic=True, crack_spacing=sth), exx, gamma, ts2_vals=(0, 0, 1))
    off = _shell_tension_then_cyclic_shear(
        lambda t: _mat(t, interlock=False), exx, gamma, ts2_vals=(0, 0, 1))
    nxy_on, nxy_off = on[-1][1], off[-1][1]
    assert nxy_on < 0.7 * nxy_off, f"interlock not load-bearing in shell: on={nxy_on} off={nxy_off}"
