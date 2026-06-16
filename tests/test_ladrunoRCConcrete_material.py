"""LadrunoRCConcrete (RC plastic-damage: ASDConcrete3D spine + MCFT compression
softening on the strength axis + Phase-2a fixed-crack aggregate-interlock shear
retention, classTag 33015) — Zone-A material battery.

Driver: a single ``stdBrick`` unit cube with 1/8-symmetry determinate restraints
develops a homogeneous *uniaxial stress* state with free lateral strain; we drive
node 2 dof 1 by displacement control and read the centroid Gauss-point stress.

The blocking A1 gate (beta on the STRENGTH axis => |sigma_c|=beta*fc' exactly,
the forbidden abscissa insertion misses it) is proven independently in
tests/_testbed/rc_shell_ref.py (numpy) and a standalone g++ build of
LadrunoRCKernel.h. Here we verify the C++ material AS INTEGRATED IN OPENSEES:

  * elastic uniaxial:                sigma_xx = E*eps_xx below cracking.
  * A2 reduce-to-ASDConcrete3D:      with beta OFF the view is byte-faithful to
    the ASDConcrete3D spine -> identical stress trajectory (tension & compression).
  * compression softening present:   uniaxial compression develops Poisson lateral
    TENSION (eps_1 = nu*|eps_axial|), so beta<1 lowers the compressive peak vs the
    beta-off run — the MCFT effect, visible end-to-end.

Plan: Ladruno_implementation/19_ladruno_rc_shell_adr.md.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}

# curved backbone shared by both materials (ASDConcrete3D -Ce/-Cs/-Cd convention;
# the kernel derives q = y/(1-d)). Compression peak fc'=30 @ 0.002 (d=0.25),
# softening to 5 @ 0.01 (d=1-5/45). Tension peak ft=3 @ 1e-4.
E, NU, KC = 30000.0, 0.2, 2.0 / 3.0
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0,   30.0,   5.0]
CD = [0.0, 0.0,    0.25,   1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010]
TS = [0.0, 3.0,    0.5]
TD = [0.0, 0.0,    1.0 - 0.5 / 5.0]


def _rc(tag, beta=False, lub_reduced=None,
        interlock=False, agg=16.0, crack_spacing=0.0, beta_sr_min=0.01):
    if lub_reduced is None:
        lub_reduced = beta
    args = ["LadrunoRCConcrete", tag, E, NU,
            "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
            "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC]
    if beta:
        args += ["-beta"]
    if lub_reduced:
        args += ["-lublinerReduced"]
    if interlock:
        args += ["-interlock", "-agg", agg, "-betaSrMin", beta_sr_min]
        if crack_spacing > 0.0:
            args += ["-crackSpacing", crack_spacing]
    ops.nDMaterial(*args)


def _asd(tag):
    ops.nDMaterial("ASDConcrete3D", tag, E, NU,
                   "-Te", *TE, "-Ts", *TS, "-Td", *TD,
                   "-Ce", *CE, "-Cs", *CS, "-Cd", *CD, "-Kc", KC)


def _build(mat_fn):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    mat_fn(1)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 1)
    ops.fix(3, 0, 0, 1)
    ops.fix(4, 1, 0, 1)
    ops.fix(5, 1, 1, 0)
    ops.fix(6, 0, 1, 0)
    ops.fix(8, 1, 0, 0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, 0.25, 0.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 100, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _run(mat_fn, eps_target, nsteps):
    """Monotonic ramp of eps_xx (node 2 dof 1) to eps_target. Returns
    [(eps_xx, sig_xx), ...]."""
    _build(mat_fn)
    out = []
    d = eps_target / nsteps
    ops.integrator("DisplacementControl", 2, 1, d)
    for _ in range(nsteps):
        assert ops.analyze(1) == 0, f"analyze failed heading to eps={eps_target}"
        ops.eleResponse(1, "forces")                      # set lazy strain
        sig = list(ops.eleResponse(1, "stresses"))[0:6]   # centroid GP
        out.append((ops.nodeDisp(2, 1), sig[0]))
    return out


# --------------------------------------------------------------------------
@pytest.mark.t0m
def test_elastic_uniaxial():
    res = _run(lambda t: _rc(t, beta=False), 1.0e-5, 5)
    eps, sig = res[-1]
    assert abs(sig - E * eps) <= 1.0e-4 * E * abs(eps) + 1.0e-8, (eps, sig, E * eps)


# --------------------------------------------------------------------------
@pytest.mark.t1
def test_reduce_to_asdconcrete3d_tension():
    """beta OFF -> LadrunoRCConcrete must track the ASDConcrete3D spine exactly,
    through the tensile peak + softening."""
    rc = _run(lambda t: _rc(t, beta=False), 5.0e-4, 60)
    asd = _run(_asd, 5.0e-4, 60)
    assert len(rc) == len(asd)
    for (e1, s1), (e2, s2) in zip(rc, asd):
        tol = 1.0e-5 * (abs(s2) if abs(s2) > 1.0 else 1.0)
        assert abs(s1 - s2) <= tol, f"eps={e1:.2e}: RC {s1} vs ASD {s2}"


@pytest.mark.t1
def test_reduce_to_asdconcrete3d_compression():
    """Same identity gate on the compression branch (past the fc' peak)."""
    rc = _run(lambda t: _rc(t, beta=False), -6.0e-3, 120)
    asd = _run(_asd, -6.0e-3, 120)
    assert len(rc) == len(asd)
    for (e1, s1), (e2, s2) in zip(rc, asd):
        tol = 1.0e-5 * (abs(s2) if abs(s2) > 1.0 else 1.0)
        assert abs(s1 - s2) <= tol, f"eps={e1:.2e}: RC {s1} vs ASD {s2}"


# --------------------------------------------------------------------------
def _homog(mat_fn, exx, eyy, ezz, nsteps=100):
    """Impose a HOMOGENEOUS strain (exx,eyy,ezz; no shear) on the unit cube by
    prescribing every nodal displacement u=(exx*x, eyy*y, ezz*z), ramped to full.
    Returns the centroid-GP stress (6-vec) at the final state."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    mat_fn(1)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 1, exx * x)
        ops.sp(t, 2, eyy * y)
        ops.sp(t, 3, ezz * z)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e14, 1.0e14)   # keep all DOFs (no 0-equation crash)
    ops.test("NormDispIncr", 1.0e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    for _ in range(nsteps):
        assert ops.analyze(1) == 0, "homogeneous-strain analyze failed"
        ops.eleResponse(1, "forces")
    return list(ops.eleResponse(1, "stresses"))[0:6]


def _min_principal(sig):
    import numpy as np
    M = np.array([[sig[0], sig[3], sig[5]],
                  [sig[3], sig[1], sig[4]],
                  [sig[5], sig[4], sig[2]]])
    return float(np.min(np.linalg.eigvalsh(M)))


@pytest.mark.t1
def test_compression_softening_biaxial():
    """The MCFT gate, end-to-end in OpenSees. Impose a biaxial state: transverse
    membrane tension eps_yy=+3e-3 (=> beta=1/(0.8+170*3e-3)=0.7634) with axial
    compression eps_xx=-2e-3, eps_zz=0. With lublinerReduced OFF, beta scales ONLY the
    assembled compressive cone and does NOT touch the damage evolution, so at the same
    end-state the most-compressive principal stress ratio (beta-ON / beta-OFF) must
    equal beta(eps_yy) exactly. (The closed-form A1 gate is also proven in
    tests/_testbed/rc_shell_ref.py + a standalone g++ build of the kernel.)"""
    e1 = 3.0e-3
    beta = 1.0 / (0.8 + 170.0 * e1)
    sig_off = _homog(lambda t: _rc(t, beta=False), -2.0e-3, e1, 0.0)
    sig_on = _homog(lambda t: _rc(t, beta=True, lub_reduced=False), -2.0e-3, e1, 0.0)
    sc_off = _min_principal(sig_off)
    sc_on = _min_principal(sig_on)
    assert sc_off < -1.0, f"no compression developed: off={sc_off}"
    ratio = sc_on / sc_off
    assert ratio == pytest.approx(beta, rel=1.0e-3), f"ratio {ratio} != beta {beta}"


# ==========================================================================
#  Phase 2a — fixed-crack aggregate-interlock shear retention
# ==========================================================================
def _path_tension_then_shear(mat_fn, exx, gamma, n1=40, n2=60, gamma0=0.0):
    """Two-stage NON-PROPORTIONAL homogeneous path on the unit cube via disjoint DOFs:
        stage 1 (pseudo-time 0->1): ramp u_x = exx*X (eps_xx) AND u_y = gamma0*X
                                    => crack freezes at the principal dir of (exx,0,gamma0)
                                       (axis-aligned when gamma0=0; oblique otherwise).
        stage 2 (pseudo-time 1->2): HOLD tension, ramp shear u_y -> gamma*X
                                    => the principal dir rotates, developing shear on the
                                       FROZEN crack plane (this is what engages interlock;
                                       a single proportional ramp would stay coaxial = no
                                       crack-plane shear to cap).
    Tension on dof-1 (∝X), shear on dof-2 (∝X, the non-symmetric simple shear u_y=gxy*X
    giving engineering gxy) -> disjoint DOFs. Returns [(gxy, sig6), ...] over stage 2."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    mat_fn(1)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Path", 1, "-time", 0.0, 1.0, 2.0, "-values", 0.0, 1.0, 1.0)  # tension: ramp, hold
    f0 = (gamma0 / gamma) if gamma != 0.0 else 0.0
    ops.timeSeries("Path", 2, "-time", 0.0, 1.0, 2.0, "-values", 0.0, f0, 1.0)   # shear: ->g0, ->gamma
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 1, exx * x)     # u_x = exx*X  -> eps_xx
        ops.sp(t, 3, 0.0)         # u_z = 0
    ops.pattern("Plain", 2, 2)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 2, gamma * x)   # u_y = gamma*X -> engineering gxy = gamma
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e14, 1.0e14)
    ops.test("NormDispIncr", 1.0e-8, 200, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / n1)
    ops.analysis("Static")
    for _ in range(n1):                                   # stage 1: pure tension
        assert ops.analyze(1) == 0, "tension stage analyze failed"
    out = []
    ops.eleResponse(1, "forces")
    out.append((0.0, list(ops.eleResponse(1, "stresses"))[0:6]))
    ops.integrator("LoadControl", 1.0 / n2)
    for _ in range(n2):                                   # stage 2: add shear
        assert ops.analyze(1) == 0, "shear stage analyze failed"
        ops.eleResponse(1, "forces")
        sig = list(ops.eleResponse(1, "stresses"))[0:6]
        gxy = gamma * (ops.getTime() - 1.0)              # ts2 factor = pseudo-time - 1
        out.append((gxy, sig))
    return out


def _vci_max(sqrt_fc, w, a_g):
    return 0.18 * sqrt_fc / (0.31 + 24.0 * w / (a_g + 16.0))


@pytest.mark.t1
def test_interlock_off_matches_asd_on_shear_path():
    """With -interlock OFF the membrane-shear path must stay byte-faithful to the
    ASDConcrete3D spine (extends the Phase-1 reduce gate to a sheared state)."""
    rc = _path_tension_then_shear(lambda t: _rc(t, beta=False, interlock=False), 5.0e-4, 2.0e-3)
    asd = _path_tension_then_shear(_asd, 5.0e-4, 2.0e-3)
    assert len(rc) == len(asd)
    for (g1, s1), (g2, s2) in zip(rc, asd):
        for k in range(6):
            tol = 1.0e-5 * (abs(s2[k]) if abs(s2[k]) > 1.0 else 1.0)
            assert abs(s1[k] - s2[k]) <= tol, f"gxy={g1:.2e} comp{k}: RC {s1[k]} vs ASD {s2[k]}"


@pytest.mark.t1
def test_interlock_shear_cap_closed_form_and_oracle():
    """-interlock ON: the membrane shear on the frozen (x-normal) crack saturates at
    the MCFT crack-shear limit v_ci,max = 0.18*sqrt(fc')/(0.31+24w/(a_g+16)),
    w = eps_n*s_theta. The crack is axis-aligned so global sigma_xy == tau_ci, which
    is backbone-independent -> we assert BOTH the closed form and the numpy oracle."""
    import numpy as np
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "_testbed"))
    import rc_shell_ref as ref

    exx, gamma, sth, a_g = 5.0e-4, 2.0e-3, 50.0, 16.0
    rc = _path_tension_then_shear(
        lambda t: _rc(t, beta=False, interlock=True, agg=a_g, crack_spacing=sth),
        exx, gamma, n1=40, n2=80)

    sqrt_fc = np.sqrt(max(CS))               # fcmax = 30
    w = exx * sth                            # crack normal along x => eps_n = eps_xx
    vci = _vci_max(sqrt_fc, w, a_g)

    # C++ shear saturates at v_ci,max at the end of the shear ramp
    tau_final = rc[-1][1][3]
    assert tau_final == pytest.approx(vci, rel=2.0e-3), f"tau_xy {tau_final} != v_ci,max {vci}"

    # ... and never overshoots the cap anywhere on the path
    for g, sig in rc:
        assert sig[3] <= vci * (1.0 + 2.0e-3), f"overshoot at gxy={g:.2e}: {sig[3]} > {vci}"

    # Independent numpy-oracle confirmation that the SAME v_ci,max cap is reached.
    # (Formulation B caps the SMEARED shear, so only the CAPPED value sigma_xy==v_ci,max
    # is backbone-independent; sub-cap sigma_xy is the damaged smeared shear and differs
    # between the raw-backbone oracle and the adjusted-backbone kernel — so we compare the
    # plateau, not the path.)
    P = ref.reference_params(beta_on=False)
    P.interlock_on = True
    P.agg_size = a_g
    P.crack_spacing = sth
    st = ref.State()
    for k in range(1, 41):                    # oracle stage 1: pure tension
        s = np.zeros(6); s[0] = exx * k / 40.0
        _, _, st = ref.compute(P, st, s)
    assert st.cracked >= 0.5 and st.crackC == pytest.approx(1.0), "oracle crack not x-aligned"
    sig_o, info_o, _ = ref.compute(P, st, np.array([exx, 0.0, 0.0, gamma, 0.0, 0.0]))
    assert sig_o[3] == pytest.approx(vci, rel=2.0e-3), f"oracle tau {sig_o[3]} != v_ci,max {vci}"
    assert info_o['beta_sr'] < 1.0, "oracle cap did not engage"


@pytest.mark.t1
def test_interlock_ablation_shear_is_load_bearing():
    """Ablation: at large post-crack shear the interlock cap must make |tau_xy| with
    interlock strictly LOWER than the smeared (interlock-off) shear -> the retention
    term is load-bearing, not cosmetic. Also asserts the crack actually formed."""
    exx, gamma, sth = 5.0e-4, 3.0e-3, 50.0
    on = _path_tension_then_shear(
        lambda t: _rc(t, beta=False, interlock=True, crack_spacing=sth), exx, gamma)
    off = _path_tension_then_shear(
        lambda t: _rc(t, beta=False, interlock=False), exx, gamma)
    tau_on, tau_off = on[-1][1][3], off[-1][1][3]
    assert tau_on < 0.5 * tau_off, f"interlock not load-bearing: on={tau_on} off={tau_off}"


@pytest.mark.t1
def test_interlock_offaxis_crack_rotation():
    """OFF-AXIS crack: stage 1 freezes the crack at an OBLIQUE angle (tension + initial
    shear gamma0), stage 2 ramps shear further so the principal direction rotates and
    develops shear ON the frozen oblique plane. This exercises the full cs!=0 rotate /
    cap / rotate-back algebra that the axis-aligned tests collapse to identity. Asserts
    (a) the reported crack normal is genuinely oblique and matches the stage-1 principal
    dir, and (b) the crack-plane shear PROJECTED from the actual stress saturates at
    v_ci,max (backbone-independent) -> rotation kernel + cap jointly correct off-axis."""
    import numpy as np
    exx, gamma0, gamma, sth, a_g = 5.0e-4, 2.0e-3, 6.0e-3, 100.0, 16.0
    rc = _path_tension_then_shear(
        lambda t: _rc(t, beta=False, interlock=True, agg=a_g, crack_spacing=sth),
        exx, gamma, n1=40, n2=80, gamma0=gamma0)
    sig = rc[-1][1]

    cs_resp = list(ops.eleResponse(1, "material", "1", "crackState"))
    assert cs_resp[0] >= 0.5, "no crack formed"
    c, s = cs_resp[1], cs_resp[2]
    c2, s2, cs = c * c, s * s, c * s

    # the frozen normal must be the stage-1 principal-strain dir (exx,0,gamma0), and oblique
    eps2 = np.array([[exx, 0.5 * gamma0], [0.5 * gamma0, 0.0]])
    w_, V_ = np.linalg.eigh(eps2)
    n = V_[:, int(np.argmax(w_))]
    assert abs(c * n[0] + s * n[1]) == pytest.approx(1.0, abs=2.0e-3), (
        f"crack normal {(c, s)} != stage-1 principal dir {tuple(n)}")
    assert 0.15 < abs(s) < 0.99, f"crack not oblique: n=({c:.3f},{s:.3f})"

    # crack-plane shear projected from the actual stress saturates at v_ci,max
    en = exx * c2 + gamma * cs                          # eyy = 0
    vci = _vci_max(np.sqrt(max(CS)), max(0.0, en) * sth, a_g)
    tau_nt = cs * (sig[1] - sig[0]) + (c2 - s2) * sig[3]
    assert abs(tau_nt) == pytest.approx(vci, rel=5.0e-3), (
        f"off-axis crack-plane shear {abs(tau_nt)} != v_ci,max {vci} (n=({c:.3f},{s:.3f}))")
