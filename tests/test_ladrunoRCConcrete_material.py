"""LadrunoRCConcrete (RC plastic-damage: ASDConcrete3D spine + MCFT compression
softening on the strength axis, classTag 33014) — Zone-A material battery.

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


def _rc(tag, beta=False, lub_reduced=None):
    if lub_reduced is None:
        lub_reduced = beta
    args = ["LadrunoRCConcrete", tag, E, NU,
            "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
            "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC]
    if beta:
        args += ["-beta"]
    if lub_reduced:
        args += ["-lublinerReduced"]
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
