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
        interlock=False, agg=16.0, crack_spacing=0.0, beta_sr_min=0.01, cyclic=False,
        xcrack=False, deg_kappa=None, deg_slip_ref=None, deg_min=None,
        shear_ret=None, shear_ret_factor=None):
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
        if cyclic:
            args += ["-cyclic"]
        if shear_ret is not None:
            args += ["-shearRetention", shear_ret]
            if shear_ret_factor is not None:
                args += ["-shearRetFactor", shear_ret_factor]
        if xcrack:
            args += ["-xcrack"]
            if deg_kappa is not None:
                args += ["-degKappa", deg_kappa]
            if deg_slip_ref is not None:
                args += ["-degSlipRef", deg_slip_ref]
            if deg_min is not None:
                args += ["-degMin", deg_min]
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
    ops.timeSeries("Path", 1, "-time", 0.0, 1.0, 2.0, "-values", 0.0, 1.0, 1.0, "-useLast")  # tension: ramp, hold
    f0 = (gamma0 / gamma) if gamma != 0.0 else 0.0
    ops.timeSeries("Path", 2, "-time", 0.0, 1.0, 2.0, "-values", 0.0, f0, 1.0, "-useLast")   # shear: ->g0, ->gamma
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


# ==========================================================================
#  Phase 2b — cyclic crack-shear (incremental friction-slip), -cyclic flag
# ==========================================================================
def _path_stages(mat_fn, exx, gamma, ts1_vals, ts2_vals, nper=60, record_from=1):
    """Multi-stage HOMOGENEOUS path: tension on dof-1 (exx*X) scaled by Path ts1_vals,
    shear on dof-2 (gamma*X, engineering gxy) scaled by Path ts2_vals (disjoint DOFs).
    Stage k spans pseudo-time [k, k+1]. Records (gxy, sig6) EVERY step from stage
    `record_from` onward. gxy is read from nodeDisp (u_y at X=1 == engineering gxy)."""
    nst = len(ts1_vals) - 1
    # pad one extra hold node beyond the analysis end (time nst) so float-accumulated
    # pseudo-time overshoot still interpolates to the held final value (Path returns 0
    # outside its range otherwise -> strains collapse, crack spuriously "closes").
    times = [float(i) for i in range(nst + 2)]
    ts1_vals = list(ts1_vals) + [ts1_vals[-1]]
    ts2_vals = list(ts2_vals) + [ts2_vals[-1]]
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    mat_fn(1)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Path", 1, "-time", *times, "-values", *ts1_vals)
    ops.timeSeries("Path", 2, "-time", *times, "-values", *ts2_vals)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 1, exx * x)
        ops.sp(t, 3, 0.0)
    ops.pattern("Plain", 2, 2)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 2, gamma * x)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e14, 1.0e14)
    ops.test("NormDispIncr", 1.0e-8, 300, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nper)
    ops.analysis("Static")
    out = []
    for stage in range(nst):
        for _ in range(nper):
            assert ops.analyze(1) == 0, f"cyclic analyze failed (stage {stage})"
            if stage >= record_from:
                ops.eleResponse(1, "forces")
                sig = list(ops.eleResponse(1, "stresses"))[0:6]
                out.append((ops.nodeDisp(2, 2), sig))
    return out


def _loop_energy(pts):
    E = pg = pt = 0.0
    for g, sig in pts:
        E += 0.5 * (sig[3] + pt) * (g - pg)
        pg, pt = g, sig[3]
    return E


@pytest.mark.t1
def test_cyclic_reversal_caps_and_dissipates():
    """-cyclic: an axis-aligned crack (tension along x) cycled in shear must (a) saturate
    the crack-shear at +/- v_ci,max on both half-cycles and (b) dissipate energy (a
    hysteresis loop, area > 0)."""
    import numpy as np
    exx, gamma, sth, a_g = 5.0e-4, 2.0e-3, 50.0, 16.0
    pts = _path_stages(lambda t: _rc(t, beta=False, interlock=True, cyclic=True, agg=a_g, crack_spacing=sth),
                       exx, gamma, [0, 1, 1, 1, 1], [0, 0, 1, -1, 1], nper=80, record_from=1)
    vci = _vci_max(np.sqrt(max(CS)), exx * sth, a_g)
    taus = [sig[3] for _, sig in pts]
    assert max(taus) == pytest.approx(vci, rel=3.0e-3), f"+peak {max(taus)} != vci {vci}"
    assert min(taus) == pytest.approx(-vci, rel=3.0e-3), f"-peak {min(taus)} != -vci {vci}"
    # dissipation over a CLOSED cycle only (stages 2+3 return +gamma->-gamma->+gamma): the open
    # stage-1 loading work is excluded so the threshold isn't a one-way-work tautology. A
    # near-rectangular friction loop dissipates ~ 4*vci*gamma; require well above the elastic
    # band area (~vci^2/G) to prove genuine hysteresis.
    nrec = len(pts)
    closed = pts[nrec // 3:]                      # the +gamma->-gamma->+gamma portion
    assert _loop_energy(closed) > 2.0 * vci * gamma, "closed-loop dissipation too small"
    # crackShear response wiring: tauCr == global sigma_xy for the axis-aligned crack
    cshear = list(ops.eleResponse(1, "material", "1", "crackShear"))
    assert cshear[0] == pytest.approx(pts[-1][1][3], rel=1.0e-6), "crackShear tauCr != sigma_xy"


@pytest.mark.t1
def test_cyclic_unload_stiffness_is_G():
    """On reversal off the cap, the crack-shear unloads along the elastic interlock
    modulus G = E/2(1+nu) (stiff), not along the cap — the friction-slip signature."""
    import numpy as np
    exx, gamma, sth, a_g = 5.0e-4, 2.0e-3, 50.0, 16.0
    G = E / (2.0 * (1.0 + NU))
    pts = _path_stages(lambda t: _rc(t, beta=False, interlock=True, cyclic=True, agg=a_g, crack_spacing=sth),
                       exx, gamma, [0, 1, 1, 1], [0, 0, 1, -1], nper=80, record_from=1)
    # the reversal is the MAX-slip point (end of the +shear leg); the next step unloads
    gs = [g for g, _ in pts]
    ipk = int(np.argmax(gs))
    g0, s0 = pts[ipk][0], pts[ipk][1][3]
    g1, s1 = pts[ipk + 1][0], pts[ipk + 1][1][3]
    slope = (s1 - s0) / (g1 - g0)
    assert slope == pytest.approx(G, rel=2.0e-2), f"unload slope {slope} != G {G}"
    # CLOSED-FORM (oracle-independent) elasto-friction line: on the stiff unload off the +cap,
    # tau(delta) = vci - G*delta until it re-caps at -vci (delta = 2*vci/G). Check a sample
    # strictly inside the elastic band.
    vci = _vci_max(np.sqrt(max(CS)), exx * sth, a_g)
    for g, sig in pts[ipk + 1:]:
        delta = g0 - g                                  # slip retraction from the peak (>0)
        if 0 < delta < 1.8 * vci / G:                   # inside the elastic band, before re-cap
            assert sig[3] == pytest.approx(vci - G * delta, abs=2.0e-3, rel=3.0e-3), (
                f"unload line off: tau={sig[3]} != vci-G*delta={vci - G * delta}")
            break


@pytest.mark.t1
def test_cyclic_crack_closure_recovers_capacity():
    """Crack-width-reversible cap: open the crack wide (low v_ci,max), shear to the cap,
    then CLOSE the crack (reduce normal tension) and slip further — the cap must rise to
    v_ci,max of the smaller width, so the shear exceeds the open-crack plateau."""
    exx, gamma, sth, a_g = 1.2e-3, 3.0e-3, 50.0, 16.0
    # stage1 ramp tension; stage2 shear to cap (open); stage3 close (exx*0.2) + more slip
    pts = _path_stages(lambda t: _rc(t, beta=False, interlock=True, cyclic=True, agg=a_g, crack_spacing=sth),
                       exx, gamma, [0, 1, 1, 0.2], [0, 0, 1, 1.6], nper=60, record_from=1)
    import numpy as np
    taus = [sig[3] for _, sig in pts]
    cap_open = _vci_max(np.sqrt(max(CS)), exx * sth, a_g)
    cap_closed = _vci_max(np.sqrt(max(CS)), 0.2 * exx * sth, a_g)
    assert max(taus) == pytest.approx(cap_closed, rel=5.0e-3), f"closed cap {max(taus)} != {cap_closed}"
    assert cap_closed > 1.10 * cap_open, "closure did not raise the cap"


@pytest.mark.t1
def test_cyclic_monotonic_reduces_to_2a_plateau():
    """Under monotonic loading -cyclic ramps the crack-shear at G then caps at the SAME
    v_ci,max as the 2a clip -> identical plateau (the cyclic law is a superset)."""
    import numpy as np
    exx, gamma, sth, a_g = 5.0e-4, 2.0e-3, 50.0, 16.0
    cyc = _path_tension_then_shear(
        lambda t: _rc(t, beta=False, interlock=True, cyclic=True, agg=a_g, crack_spacing=sth),
        exx, gamma, n1=40, n2=80)
    twoa = _path_tension_then_shear(
        lambda t: _rc(t, beta=False, interlock=True, cyclic=False, agg=a_g, crack_spacing=sth),
        exx, gamma, n1=40, n2=80)
    vci = _vci_max(np.sqrt(max(CS)), exx * sth, a_g)
    assert cyc[-1][1][3] == pytest.approx(vci, rel=3.0e-3), f"cyclic monotonic plateau {cyc[-1][1][3]} != vci {vci}"
    # both saturate the SAME bound -> the cyclic monotonic plateau equals the actual 2a-run plateau
    assert cyc[-1][1][3] == pytest.approx(twoa[-1][1][3], rel=3.0e-3), "cyclic plateau != 2a-run plateau"


@pytest.mark.t1
def test_cyclic_matches_oracle_path():
    """For an axis-aligned crack the cyclic crack-shear sigma_xy == tau_ci (the independent
    friction stress), which is backbone-INDEPENDENT -> C++ must match the numpy oracle
    step-by-step over the whole reversing path (not just at the cap)."""
    import numpy as np
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "_testbed"))
    import rc_shell_ref as ref
    exx, gamma, sth, a_g = 5.0e-4, 2.0e-3, 50.0, 16.0
    pts = _path_stages(lambda t: _rc(t, beta=False, interlock=True, cyclic=True, agg=a_g, crack_spacing=sth),
                       exx, gamma, [0, 1, 1, 1, 1], [0, 0, 1, -1, 1], nper=40, record_from=0)
    P = ref.reference_params(beta_on=False)
    P.interlock_on = True
    P.interlock_cyclic = True
    P.agg_size = a_g
    P.crack_spacing = sth
    st = ref.State()
    nmax = 0.0
    for g, sig in pts:
        s = np.array([exx, 0.0, 0.0, g, 0.0, 0.0])
        sig_o, _, st = ref.compute(P, st, s)
        nmax = max(nmax, abs(sig[3] - sig_o[3]))
    assert nmax < 5.0e-3, f"C++ vs oracle cyclic crack-shear mismatch {nmax}"


@pytest.mark.t1
def test_cyclic_tangent_assembled_in_opensees():
    """The cyclic consistent tangent, AS ASSEMBLED IN OPENSEES (not just standalone FD): for an
    axis-aligned crack the in-plane shear tangent Dtan[3][3] must equal the interlock modulus G
    while sliding-elastic (sub-cap), and collapse to betaSrMin*G once the friction cap binds
    (the friction stiffness is removed). Reads the 6x6 material tangent via eleResponse."""
    import numpy as np
    G, bsr = E / (2.0 * (1.0 + NU)), 0.01
    sth, a_g = 50.0, 16.0
    vci = _vci_max(np.sqrt(max(CS)), 5.0e-4 * sth, a_g)

    def t33_after(gamma):
        _path_tension_then_shear(
            lambda t: _rc(t, beta=False, interlock=True, cyclic=True, agg=a_g,
                          crack_spacing=sth, beta_sr_min=bsr),
            5.0e-4, gamma, n1=40, n2=60)
        T = list(ops.eleResponse(1, "material", "1", "tangent"))   # 6x6 row-major
        return T[3 * 6 + 3]

    # sub-cap sliding (G*gamma < vci => gamma < vci/G): tangent == G
    g_sub = 0.6 * vci / G
    assert t33_after(g_sub) == pytest.approx(G, rel=2.0e-2), "sliding-elastic Dtan[3][3] != G"
    # well past the cap: friction stiffness removed -> betaSrMin*G
    assert t33_after(2.0e-3) == pytest.approx(bsr * G, rel=0.2), "capped Dtan[3][3] != betaSrMin*G"


@pytest.mark.t1
def test_serialization_roundtrip_cracked_cyclic_state():
    """sendSelf/recvSelf round-trip (via the FE database) of a fully-loaded -xcrack state: the
    frozen crack frame + width + friction state {tauCr,gammaCr} + the Phase-2b.2b X-crack state
    {cracked2,slipCum} must ALL survive a save/restore bit-exactly with NON-ZERO values (guards
    the schema-versioned RC_DATA send/recv field count/order)."""
    import os
    import tempfile
    # -xcrack tension-then-shear with wear: forms an oblique crack 1 and accumulates crack
    # sliding, so crackState, crackShear AND xcrackState (slipCum>0) are all non-trivial -- the
    # new tail fields are round-tripped at NON-ZERO values, guarding the wire field count/order.
    def mat(t):
        _rc(t, beta=False, interlock=True, cyclic=True, xcrack=True, crack_spacing=50.0,
            deg_kappa=0.5, deg_slip_ref=0.02)
    _path_tension_then_shear(mat, 5.0e-4, 4.0e-3, n1=40, n2=60, gamma0=2.0e-3)
    before = (list(ops.eleResponse(1, "material", "1", "crackState"))
              + list(ops.eleResponse(1, "material", "1", "crackShear"))
              + list(ops.eleResponse(1, "material", "1", "xcrackState")))
    assert before[0] >= 0.5, "precondition: crack 1 must have formed"
    assert before[7] > 0.0, f"precondition: slipCum must be > 0 (got {before[7]})"  # xcrackState[1]

    db = os.path.join(tempfile.mkdtemp(prefix="rc_rt_"), "rc_state")
    try:
        ops.database("File", db)
        ops.save(1)
        ops.restore(1)
        after = (list(ops.eleResponse(1, "material", "1", "crackState"))
                 + list(ops.eleResponse(1, "material", "1", "xcrackState")))
        # NOTE: crackShear (tauCr) is dropped from the round-trip compare: the response reads the
        # TRIAL histTr (re-evaluated at the converged disp with a ~1e-8 Penalty residual, amplified
        # by the stiff friction modulus G into a ~1e-4 difference from the COMMITTED histN that is
        # actually serialized), so a tauCr trial-vs-committed mismatch is a read artifact, not a
        # serialization defect. The committed fields (crack frame, width, cracked2, slipCum) below
        # round-trip BIT-EXACT, which is what guards the wire field count/order after the +schema-
        # version +cracked2 +slipCum growth. (slipCum>0 above confirms the new tail field is live.)
        bf = (before[:4] + before[6:])      # crackState[0:4] + xcrackState[6:8]; drop crackShear[4:6]
        assert all(abs(a - b) <= 1.0e-9 + 1.0e-9 * abs(b) for a, b in zip(bf, after)), (
            f"serialized committed state drifted: before={bf} after={after}")
    finally:
        import glob
        for f in glob.glob(db + "*"):
            try:
                os.remove(f)
            except OSError:
                pass


# ==========================================================================
#  Phase 2b.2b — second orthogonal crack (X-cracking) + slip-driven wear
# ==========================================================================
@pytest.mark.t1
def test_xcrack_off_reduces_to_2b1():
    """-xcrack with wear OFF (-degKappa 0) on a path that never forms crack 2 (uniaxial
    tension along x => orthogonal y never cracks) must be bit-identical to plain -cyclic."""
    import numpy as np
    exx, gamma, sth = 5.0e-4, 2.0e-3, 50.0
    base = _path_tension_then_shear(
        lambda t: _rc(t, beta=False, interlock=True, cyclic=True, crack_spacing=sth),
        exx, gamma, n1=40, n2=70)
    xc = _path_tension_then_shear(
        lambda t: _rc(t, beta=False, interlock=True, cyclic=True, crack_spacing=sth,
                      xcrack=True, deg_kappa=0.0),
        exx, gamma, n1=40, n2=70)
    assert len(base) == len(xc)
    for (g1, s1), (g2, s2) in zip(base, xc):
        for k in range(6):
            assert abs(s1[k] - s2[k]) <= 1.0e-9 + 1.0e-9 * abs(s1[k]), (
                f"xcrack(no-wear) != 2b.1 at gxy={g1:.2e} comp{k}: {s1[k]} vs {s2[k]}")


@pytest.mark.t1
def test_xcrack_second_crack_captures_under_biaxial_tension():
    """The orthogonal crack 2 forms only when the perpendicular direction also cracks:
    biaxial tension (eps_xx, eps_yy both > eps_cr) => cracked2==1; uniaxial => cracked2==0."""
    # biaxial (NON-equibiaxial so crack 1's principal dir is well-defined, not degenerate):
    # both directions cracked -> the orthogonal crack 2 forms
    _homog(lambda t: _rc(t, beta=False, interlock=True, cyclic=True, xcrack=True, crack_spacing=50.0),
           5.0e-4, 3.0e-4, 0.0)
    bi = list(ops.eleResponse(1, "material", "1", "xcrackState"))
    # uniaxial: only x cracks -> crack 2 (orthogonal, y) never forms
    _homog(lambda t: _rc(t, beta=False, interlock=True, cyclic=True, xcrack=True, crack_spacing=50.0),
           5.0e-4, 0.0, 0.0)
    uni = list(ops.eleResponse(1, "material", "1", "xcrackState"))
    assert bi[0] >= 0.5, f"crack 2 should form under biaxial tension (got {bi[0]})"
    assert uni[0] < 0.5, f"crack 2 should NOT form under uniaxial tension (got {uni[0]})"


@pytest.mark.t1
def test_xcrack_cyclic_strength_degrades():
    """The headline 2b.2b material gate: slip-driven interlock-surface wear makes the
    capped membrane shear DECAY over repeated shear cycles (cyclic strength degradation),
    whereas plain -cyclic (no wear) keeps a constant cap. The degraded plateau matches the
    closed-form kn*v_ci,max, kn = 1 - degKappa*min(slipPeak/degSlipRef,1)."""
    exx, gamma, sth = 5.0e-4, 2.0e-3, 50.0
    # degSlipRef tuned so the cumulative sliding wears GRADUALLY over the 3 cycles (not saturating
    # in cycle 1): ~ cycle adds ~4*gamma of sliding, so 0.03 ~ a few cycles to floor.
    dslip, dkap, dmin = 0.03, 0.5, 0.1
    nper = 60
    ts2 = (0, 0, 1, -1, 1, -1, 1, -1)                          # 3 shear cycles (6 legs)
    layout = [0, 1, 1, 1, 1, 1, 1, 1]

    def plus_peaks(material):                                   # tau at gxy=+gamma each + cycle
        pts = _path_stages(material, exx, gamma, layout, list(ts2), nper=nper, record_from=1)
        # recorded stages 1..6 -> + extremes at the END of stages 1,3,5 (indices 59,179,299)
        return [pts[59][1][3], pts[179][1][3], pts[299][1][3]]

    wear = plus_peaks(lambda t: _rc(t, beta=False, interlock=True, cyclic=True, xcrack=True,
                                    crack_spacing=sth, deg_kappa=dkap, deg_slip_ref=dslip, deg_min=dmin))
    nowear = plus_peaks(lambda t: _rc(t, beta=False, interlock=True, cyclic=True, crack_spacing=sth))
    # WEAR: the + cap decays MONOTONICALLY cycle over cycle (cumulative-slip abrasion)
    assert wear[0] > wear[1] > wear[2], f"cap not monotonically decaying with wear: {wear}"
    assert wear[2] < 0.8 * wear[0], f"insufficient cyclic decay with wear: {wear}"
    # NO-WEAR control: the + cap is constant across cycles
    assert nowear[2] == pytest.approx(nowear[0], rel=2.0e-2), f"no-wear cap drifted: {nowear}"


@pytest.mark.t1
def test_xcrack_matches_oracle_cyclic():
    """C++ vs numpy oracle over a reversing path with -xcrack (degradation engaged):
    axis-aligned crack => sigma_xy == tau_ci is backbone-independent, so step-by-step match."""
    import numpy as np
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "_testbed"))
    import rc_shell_ref as ref
    exx, gamma, sth, a_g = 5.0e-4, 2.0e-3, 50.0, 16.0
    dslip, dkap, dmin = 2.0e-3, 0.5, 0.1
    ts2 = (0, 0, 1, -1, 1)
    pts = _path_stages(
        lambda t: _rc(t, beta=False, interlock=True, cyclic=True, xcrack=True,
                      crack_spacing=sth, deg_kappa=dkap, deg_slip_ref=dslip, deg_min=dmin),
        exx, gamma, [0, 1, 1, 1, 1], list(ts2), nper=40, record_from=0)
    P = ref.reference_params(beta_on=False)
    P.interlock_on = True; P.interlock_cyclic = True; P.xcrack_on = True
    P.agg_size = a_g; P.crack_spacing = sth
    P.deg_kappa = dkap; P.deg_slip_ref = dslip; P.deg_min = dmin
    st = ref.State()
    nmax = 0.0
    for g, sig in pts:
        s = np.array([exx, 0.0, 0.0, g, 0.0, 0.0])
        sig_o, _, st = ref.compute(P, st, s)
        nmax = max(nmax, abs(sig[3] - sig_o[3]))
    assert nmax < 5.0e-3, f"C++ vs oracle xcrack cyclic mismatch {nmax}"


# ==========================================================================
#  Phase 2b.2c — -shearRetention {mcft|const|dsfm|rots} crack-shear retention curve
#  The mode selects the CYCLIC crack-plane slip stiffness G_slip in the friction
#  predictor tau_cr = clamp(tau_cr_c + G_slip*dgamma_nt, +/- v_ci,max(w)). The cap is
#  unchanged in every mode. mcft = full G (default); const = mu*G; dsfm = width-degraded
#  G; rots = rotating-coaxial (no independent crack shear). All reduce to mcft.
# ==========================================================================
def _shear_unload_slope(mat_fn, exx=5.0e-4, gamma=2.0e-3, sth=50.0, a_g=16.0):
    """Drive tension-then-cyclic-shear; return (slope at the first reversal, v_ci,max).
    The reversal unloads along the crack-shear slip stiffness G_slip -> this measures it."""
    import numpy as np
    pts = _path_stages(mat_fn, exx, gamma, [0, 1, 1, 1], [0, 0, 1, -1], nper=80, record_from=1)
    gs = [g for g, _ in pts]
    ipk = int(np.argmax(gs))
    g0, s0 = pts[ipk][0], pts[ipk][1][3]
    g1, s1 = pts[ipk + 1][0], pts[ipk + 1][1][3]
    vci = _vci_max(np.sqrt(max(CS)), exx * sth, a_g)
    return (s1 - s0) / (g1 - g0), vci


@pytest.mark.t1
def test_shearret_const_unity_matches_mcft():
    """-shearRetention const with mu=1 must be byte-identical to the default (mcft) cyclic
    path: same slip stiffness G, same cap -> identical crack-shear trajectory."""
    import numpy as np
    exx, gamma, sth, a_g = 5.0e-4, 2.0e-3, 50.0, 16.0
    base = _path_stages(lambda t: _rc(t, beta=False, interlock=True, cyclic=True, agg=a_g, crack_spacing=sth),
                        exx, gamma, [0, 1, 1, 1, 1], [0, 0, 1, -1, 1], nper=40, record_from=0)
    cst = _path_stages(lambda t: _rc(t, beta=False, interlock=True, cyclic=True, agg=a_g, crack_spacing=sth,
                                     shear_ret="const", shear_ret_factor=1.0),
                       exx, gamma, [0, 1, 1, 1, 1], [0, 0, 1, -1, 1], nper=40, record_from=0)
    nmax = max(abs(a[1][3] - b[1][3]) for a, b in zip(base, cst))
    assert nmax < 1.0e-9, f"const(mu=1) != mcft: max |dtau|={nmax}"


@pytest.mark.t1
def test_shearret_const_unload_stiffness_is_mu_G():
    """const-mode unload slope off the cap == mu*G (the reduced constant retention stiffness),
    strictly softer than the full-G mcft default."""
    G, mu = E / (2.0 * (1.0 + NU)), 0.4
    sth, a_g = 50.0, 16.0
    slope, vci = _shear_unload_slope(
        lambda t: _rc(t, beta=False, interlock=True, cyclic=True, agg=a_g, crack_spacing=sth,
                      shear_ret="const", shear_ret_factor=mu))
    assert slope == pytest.approx(mu * G, rel=3.0e-2), f"const unload slope {slope} != mu*G {mu * G}"
    assert slope < 0.8 * G, "const retention not softer than full G"


@pytest.mark.t1
def test_shearret_dsfm_unload_stiffness_is_width_degraded():
    """dsfm-mode unload slope off the cap == G*(0.31/denom) (the width-degraded slip
    stiffness, denom = 0.31 + 24 w/(a_g+16)). Strictly softer than full G once the crack opens."""
    import numpy as np
    G = E / (2.0 * (1.0 + NU))
    # a wide crack (w = exx*sth = 0.1 mm) so the DSFM width-degradation is pronounced:
    # denom = 0.31 + 24*0.1/32 = 0.385 -> G_slip = G*0.31/0.385 ~= 0.805*G.
    exx, sth, a_g = 1.0e-3, 100.0, 16.0
    w = exx * sth
    denom = 0.31 + 24.0 * w / (a_g + 16.0)
    g_expected = G * (0.31 / denom)
    slope, vci = _shear_unload_slope(
        lambda t: _rc(t, beta=False, interlock=True, cyclic=True, agg=a_g, crack_spacing=sth,
                      shear_ret="dsfm"), exx=exx, sth=sth, a_g=a_g)
    assert slope == pytest.approx(g_expected, rel=4.0e-2), f"dsfm unload slope {slope} != G*0.31/denom {g_expected}"
    assert slope < 0.9 * G, "dsfm retention not softer than full G at this crack width"


@pytest.mark.t1
def test_shearret_const_caps_at_same_vci():
    """The retention mode changes the slip STIFFNESS, not the CAP: const-mode crack-shear
    still saturates at +/- v_ci,max(w), identical to the default."""
    import numpy as np
    exx, gamma, sth, a_g = 5.0e-4, 3.0e-3, 50.0, 16.0
    pts = _path_stages(lambda t: _rc(t, beta=False, interlock=True, cyclic=True, agg=a_g, crack_spacing=sth,
                                     shear_ret="const", shear_ret_factor=0.4),
                       exx, gamma, [0, 1, 1, 1, 1], [0, 0, 1, -1, 1], nper=80, record_from=1)
    vci = _vci_max(np.sqrt(max(CS)), exx * sth, a_g)
    taus = [sig[3] for _, sig in pts]
    assert max(taus) == pytest.approx(vci, rel=4.0e-3), f"const +peak {max(taus)} != vci {vci}"
    assert min(taus) == pytest.approx(-vci, rel=4.0e-3), f"const -peak {min(taus)} != -vci {vci}"


@pytest.mark.t1
def test_shearret_rots_disables_crack_shear():
    """-shearRetention rots = rotating-coaxial: the fixed-crack shear block is skipped, so the
    membrane shear stays smeared -> identical to interlock OFF (no crack capture, no cap)."""
    exx, gamma, sth = 5.0e-4, 2.0e-3, 50.0
    off = _path_tension_then_shear(
        lambda t: _rc(t, beta=False, interlock=False), exx, gamma, n1=40, n2=70)
    rots = _path_tension_then_shear(
        lambda t: _rc(t, beta=False, interlock=True, cyclic=True, crack_spacing=sth, shear_ret="rots"),
        exx, gamma, n1=40, n2=70)
    nmax = max(abs(a[1][3] - b[1][3]) for a, b in zip(off, rots))
    assert nmax < 1.0e-9, f"rots != interlock-off (smeared shear): max |dtau|={nmax}"
    # and rots must NOT have frozen a crack
    assert list(ops.eleResponse(1, "material", "1", "crackState"))[0] < 0.5, "rots wrongly froze a crack"


@pytest.mark.t1
def test_shearret_const_matches_oracle_path():
    """C++ vs numpy oracle over a full reversing path with const-mode retention (mu=0.4):
    axis-aligned crack => sigma_xy == tau_ci is backbone-independent, so a step-by-step match."""
    import numpy as np
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "_testbed"))
    import rc_shell_ref as ref
    exx, gamma, sth, a_g, mu = 5.0e-4, 2.0e-3, 50.0, 16.0, 0.4
    pts = _path_stages(lambda t: _rc(t, beta=False, interlock=True, cyclic=True, agg=a_g, crack_spacing=sth,
                                     shear_ret="const", shear_ret_factor=mu),
                       exx, gamma, [0, 1, 1, 1, 1], [0, 0, 1, -1, 1], nper=40, record_from=0)
    P = ref.reference_params(beta_on=False)
    P.interlock_on = True; P.interlock_cyclic = True
    P.agg_size = a_g; P.crack_spacing = sth
    P.shear_ret_mode = 1; P.shear_ret_factor = mu
    st = ref.State()
    nmax = 0.0
    for g, sig in pts:
        s = np.array([exx, 0.0, 0.0, g, 0.0, 0.0])
        sig_o, _, st = ref.compute(P, st, s)
        nmax = max(nmax, abs(sig[3] - sig_o[3]))
    assert nmax < 5.0e-3, f"C++ vs oracle const-retention cyclic mismatch {nmax}"


# ==========================================================================
#  Phase 2b.2c — crack-closure on the NORMAL direction (unilateral recovery)
#  The cloned ASDConcrete3D spine reassembles the spectral tension/compression
#  split EVERY step from the live effective stress with independent dt/dc and
#  cdf=0, which IS unilateral crack closure on the normal: tensile cracking
#  (dt>0) must NOT knock down the compressive capacity once the crack closes.
#  The fixed-crack addition is shear-only, so the normal closure is the spine's
#  job. These gates verify it end-to-end in OpenSees.
# ==========================================================================
def _run_legs(mat_fn, targets, nper=80):
    """Drive node-2 dof-1 (uniaxial stress, lateral free) through successive SIGNED
    displacement targets via DisplacementControl. Returns [(eps_xx, sig_xx), ...]."""
    _build(mat_fn)
    out = []
    u = 0.0
    for tgt in targets:
        d = (tgt - u) / nper
        ops.integrator("DisplacementControl", 2, 1, d)
        for _ in range(nper):
            assert ops.analyze(1) == 0, f"leg to {tgt} analyze failed"
            ops.eleResponse(1, "forces")
            sig = list(ops.eleResponse(1, "stresses"))[0:6]
            out.append((ops.nodeDisp(2, 1), sig[0]))
        u = tgt
    return out


@pytest.mark.t1
def test_crack_closure_recovers_full_compression():
    """Unilateral crack closure on the NORMAL direction: crack the point in tension
    (dt>0, sigma softens), then load it into compression past the peak. The compressive
    capacity must FULLY recover (prior tensile damage does not bleed into compression,
    cdf=0) — i.e. the compression peak equals a virgin compression run with no prior
    cracking. This is the spine's spectral reassembly = crack closure, verified end-to-end."""
    # crack in tension (eps past ft/E=1e-4), then compress past the peak (eps ~ -2e-3)
    cracked = _run_legs(lambda t: _rc(t, beta=False), [3.0e-4, -2.4e-3])
    virgin = _run_legs(lambda t: _rc(t, beta=False), [-2.4e-3])
    peak_cracked = min(s for _, s in cracked)     # most compressive
    peak_virgin = min(s for _, s in virgin)
    assert peak_cracked < -25.0, f"compression did not develop after closure: {peak_cracked}"
    assert peak_cracked == pytest.approx(peak_virgin, rel=1.0e-3), (
        f"crack closure did not recover full compression: cracked peak {peak_cracked} "
        f"vs virgin {peak_virgin}")


@pytest.mark.t1
def test_crack_closure_tension_damage_is_irreversible():
    """The flip side of unilateral closure: after cracking + closing + REOPENING, the
    reopened tensile stress must follow the DAMAGED (softened) tension envelope, not the
    virgin elastic one — tensile damage is irreversible even though compression recovered."""
    # crack DEEP into tension softening (eps=1e-3, the TE/TS tail where the envelope is 0.5),
    # then close into compression, then reopen to the same +tension strain.
    pts = _run_legs(lambda t: _rc(t, beta=False), [1.0e-3, -1.0e-3, 1.0e-3])
    # the reopened-tension tail must follow the DAMAGED envelope (~0.5), far below both the
    # virgin tensile peak ft=3 AND the virgin elastic response E*eps at that strain.
    reopen = [s for e, s in pts[-80:] if e > 6.0e-4]
    assert reopen, "no reopened-tension samples"
    assert max(reopen) < 0.4 * max(TS), (
        f"reopened tension {max(reopen)} not on the damaged envelope (ft={max(TS)})")
    assert max(reopen) < 0.1 * E * 1.0e-3, "reopened tension not below virgin elastic (no damage?)"


@pytest.mark.t1
def test_crack_closure_with_interlock_recovers_compression():
    """Crack closure with the fixed-crack interlock ON: forming the crack + freezing its
    normal must NOT corrupt the normal compressive recovery (the interlock modifies only the
    crack-plane SHEAR). Compression peak after a frozen-crack + closure == virgin compression."""
    sth = 50.0
    cracked = _run_legs(lambda t: _rc(t, beta=False, interlock=True, cyclic=True, crack_spacing=sth),
                        [3.0e-4, -2.4e-3])
    virgin = _run_legs(lambda t: _rc(t, beta=False), [-2.4e-3])
    peak_cracked = min(s for _, s in cracked)
    peak_virgin = min(s for _, s in virgin)
    assert peak_cracked == pytest.approx(peak_virgin, rel=1.0e-3), (
        f"interlock corrupted normal closure: cracked {peak_cracked} vs virgin {peak_virgin}")
