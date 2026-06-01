"""LadrunoBrick (Ladruno unified 8-node hexahedron, classTag 33002) — Zone-A
element battery.

v1 ships the small-strain ``std`` + ``bbar`` formulations. The headline gate is
**bit-for-bit overlap with the upstream elements**:

  * ``-formulation std``  must reproduce ``stdBrick``   (Ed Love's Brick)
  * ``-formulation bbar`` must reproduce ``bbarBrick``  (Ed Love's BbarBrick)

to ~1e-9 on a distorted hex under mixed loading. That equivalence is our
correctness anchor (only possible because v1 is geometrically linear). On top of
that:

  * rank / zero-energy: a free element has exactly 6 rigid-body modes (18 with
    energy) for both std and bbar — no spurious mechanisms;
  * volumetric-locking relief: near-incompressible (nu=0.4999), bbar is strictly
    more flexible than std (std locks);
  * the reserved formulations (uri/eas) are refused at the factory in v1.

Plan: Ladruno_implementation/09_ladruno_brick.md.
"""
import math

import pytest

from _testbed import ops
from _testbed.fem_checks import zero_energy_mode_count, n_rigid_body_modes

pytestmark = [pytest.mark.zone_a]

# A deliberately distorted (but positive-Jacobian) hex. Standard Brick node
# order: nodes 1-4 are the z- face (CCW), nodes 5-8 the z+ face (CCW).
_NODES = {
    1: (0.00, 0.00, 0.00),
    2: (1.00, 0.10, 0.00),
    3: (1.10, 1.00, 0.10),
    4: (0.05, 0.95, 0.00),
    5: (0.00, 0.05, 1.00),
    6: (1.00, 0.00, 1.05),
    7: (1.05, 1.00, 1.10),
    8: (0.00, 1.00, 0.95),
}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]
_BASE = [1, 2, 3, 4]          # fully fixed
_TOP = [5, 6, 7, 8]           # loaded
_LOADS = {                    # mixed so the full 24x24 K is exercised
    5: (10.0, 0.0, 0.0),
    6: (0.0, 8.0, 0.0),
    7: (0.0, 0.0, -12.0),
    8: (5.0, 5.0, 5.0),
}


def _build_common(E, nu):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)


def _solve_hex(ele_name, extra_args, E=1000.0, nu=0.3):
    """Build a single hex of the given element, static-solve under _LOADS, and
    return (disps[24], stresses[48], forces[24])."""
    _build_common(E, nu)
    for n in _BASE:
        ops.fix(n, 1, 1, 1)

    ops.element(ele_name, 1, *_CONN, 1, *extra_args)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (fx, fy, fz) in _LOADS.items():
        ops.load(n, fx, fy, fz)

    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    disps = []
    for n in _CONN:
        disps.extend(ops.nodeDisp(n))
    # Read "forces" BEFORE "stresses": upstream bbarBrick has no update() and
    # sets the material trial strain only inside formResidAndTangent, so a bare
    # "stresses" read reflects the predictor (u=0) state. getResistingForce
    # ("forces") triggers that strain evaluation at the committed displacement,
    # after which "stresses" reflects the solved state for lazy-strain elements
    # too. (LadrunoBrick/Brick implement update(), so order is immaterial for
    # them.) See LEDGER_quirks: bbarBrick lazy-strain readback.
    forces = list(ops.eleResponse(1, "forces"))
    stresses = list(ops.eleResponse(1, "stresses"))
    return disps, stresses, forces


def _assert_close(a, b, label, rtol=1e-9, atol=1e-9):
    assert len(a) == len(b), f"{label}: length mismatch {len(a)} vs {len(b)}"
    for i, (x, y) in enumerate(zip(a, b)):
        tol = atol + rtol * max(abs(x), abs(y))
        assert abs(x - y) <= tol, (
            f"{label}[{i}]: {x!r} vs {y!r} (|d|={abs(x-y):.3e} > {tol:.3e})"
        )


# --------------------------------------------------------------------------
# headline regression gates
# --------------------------------------------------------------------------
def test_std_matches_stdBrick():
    """-formulation std reproduces upstream stdBrick to ~1e-9."""
    ref = _solve_hex("stdBrick", [])
    ours = _solve_hex("LadrunoBrick", ["-formulation", "std"])
    _assert_close(ours[0], ref[0], "disp")
    _assert_close(ours[1], ref[1], "stress")
    _assert_close(ours[2], ref[2], "force")


def test_bbar_matches_bbarBrick():
    """-formulation bbar reproduces upstream bbarBrick to ~1e-9."""
    ref = _solve_hex("bbarBrick", [])
    ours = _solve_hex("LadrunoBrick", ["-formulation", "bbar"])
    _assert_close(ours[0], ref[0], "disp")
    _assert_close(ours[1], ref[1], "stress")
    _assert_close(ours[2], ref[2], "force")


def test_default_formulation_is_std():
    """No -formulation flag defaults to std."""
    ref = _solve_hex("LadrunoBrick", ["-formulation", "std"])
    bare = _solve_hex("LadrunoBrick", [])
    _assert_close(bare[0], ref[0], "disp")


# --------------------------------------------------------------------------
# rank / zero-energy (no spurious mechanisms for std/bbar)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("formulation", ["std", "bbar", "uri", "eas"])
def test_free_element_has_six_rigid_body_modes(formulation):
    # For uri this is also the hourglass-stabilization gate: 1-pt integration
    # alone leaves 6 rigid + 12 hourglass = 18 zero-energy modes; correct
    # Flanagan-Belytschko hourglass control restores stiffness to the 12 so
    # only the 6 rigid-body modes remain at zero energy.
    _build_common(E=1000.0, nu=0.3)
    ops.element("LadrunoBrick", 1, *_CONN, 1, "-formulation", formulation)
    count, eigs = zero_energy_mode_count(_CONN, ndf=3)
    assert count == n_rigid_body_modes(3), (
        f"{formulation}: {count} zero-energy modes (want 6); eigs={eigs}"
    )


def test_uri_relieves_volumetric_locking():
    """1-pt reduced integration removes the over-constraint of incompressibility,
    so uri (like bbar) is more flexible than std near nu=0.5 — visible even in a
    single element (volumetric locking is a constant-mode effect)."""
    nu = 0.4999
    std = _solve_hex("LadrunoBrick", ["-formulation", "std"], nu=nu)
    uri = _solve_hex("LadrunoBrick", ["-formulation", "uri"], nu=nu)
    std_mag = math.sqrt(sum(u * u for u in std[0]))
    uri_mag = math.sqrt(sum(u * u for u in uri[0]))
    assert uri_mag > std_mag, (
        f"uri (|u|={uri_mag:.3e}) should be more flexible than std "
        f"(|u|={std_mag:.3e}) at nu={nu}"
    )


@pytest.mark.parametrize("formulation", ["std", "bbar", "uri", "eas"])
def test_uniaxial_stress_patch(formulation):
    """Patch test: a unit cube under uniform x-tractions with statically-
    determinate (6 rigid-mode) restraints develops uniform uniaxial stress
    sigma_xx = sig0 (others zero). Every formulation must reproduce it exactly.
    For uri this confirms the FB hourglass control adds exactly zero to a
    constant-strain state (gamma is orthogonal to the linear field)."""
    E, nu, sig0 = 1000.0, 0.3, 2.0
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    cube = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
            5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
    for t, c in cube.items():
        ops.node(t, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    # 1/8-symmetry restraints, all consistent with the exact field
    # u=((sig0/E)x, -nu(sig0/E)y, -nu(sig0/E)z): each min-face pinned in its normal,
    # so x=0 reacts the tension and lateral contraction stays free.
    ops.fix(1, 1, 1, 1)   # (0,0,0): x0,y0,z0
    ops.fix(2, 0, 1, 1)   # (1,0,0): y0,z0
    ops.fix(3, 0, 0, 1)   # (1,1,0): z0
    ops.fix(4, 1, 0, 1)   # (0,1,0): x0,z0
    ops.fix(5, 1, 1, 0)   # (0,0,1): x0,y0
    ops.fix(6, 0, 1, 0)   # (1,0,1): y0
    ops.fix(8, 1, 0, 0)   # (0,1,1): x0
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1, "-formulation", formulation)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):            # uniform traction sig0 on the unit x=1 face
        ops.load(n, sig0 / 4.0, 0.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    ops.eleResponse(1, "forces")                      # ensure strain is set (lazy elements)
    sig = list(ops.eleResponse(1, "stresses"))[0:6]   # centroid (Gauss point 1)
    sig_exact = [sig0, 0.0, 0.0, 0.0, 0.0, 0.0]
    _assert_close(sig, sig_exact, f"patch[{formulation}]", rtol=1e-7, atol=1e-7)


# --------------------------------------------------------------------------
# volumetric-locking relief: bbar is strictly more flexible near-incompressible
# --------------------------------------------------------------------------
def test_bbar_relieves_volumetric_locking():
    nu = 0.4999
    std = _solve_hex("LadrunoBrick", ["-formulation", "std"], nu=nu)
    bbar = _solve_hex("LadrunoBrick", ["-formulation", "bbar"], nu=nu)
    std_mag = math.sqrt(sum(u * u for u in std[0]))
    bbar_mag = math.sqrt(sum(u * u for u in bbar[0]))
    # std locks (over-stiff) -> smaller displacement magnitude than bbar.
    assert bbar_mag > std_mag, (
        f"bbar (|u|={bbar_mag:.3e}) should be more flexible than std "
        f"(|u|={std_mag:.3e}) at nu={nu}"
    )


# --------------------------------------------------------------------------
# reserved formulations / hourglass flavours refused in this build
# --------------------------------------------------------------------------
def _refused(*extra):
    _build_common(E=1000.0, nu=0.3)
    try:
        ops.element("LadrunoBrick", 1, *_CONN, 1, *extra)
    except Exception:
        pass  # raising is an acceptable refusal
    tags = ops.getEleTags() or []
    if isinstance(tags, int):
        tags = [tags]
    return 1 not in tags


def test_eas_relieves_volumetric_locking():
    """eas (bbar + statically-condensed enhanced strain) is the general-nu
    element: like bbar/uri it removes the incompressibility over-constraint, so
    it is strictly more flexible than std near nu=0.5. (Bending/all-nu accuracy
    vs SSPbrick is gated in test_ladrunoBrick_bending.py.)"""
    nu = 0.4999
    std = _solve_hex("LadrunoBrick", ["-formulation", "std"], nu=nu)
    eas = _solve_hex("LadrunoBrick", ["-formulation", "eas"], nu=nu)
    std_mag = math.sqrt(sum(u * u for u in std[0]))
    eas_mag = math.sqrt(sum(u * u for u in eas[0]))
    assert eas_mag > std_mag, (
        f"eas (|u|={eas_mag:.3e}) should be more flexible than std "
        f"(|u|={std_mag:.3e}) at nu={nu}"
    )


# --------------------------------------------------------------------------
# uri -hourglass viscous: rate-form damping, explicit-only. Patch/rank cannot
# see it (it is velocity-dependent and adds no stiffness); an explicit dynamic
# run can — it must accept the option and stay stable (no hourglass blow-up).
# --------------------------------------------------------------------------
def test_uri_viscous_hourglass_accepted():
    """The factory now builds -hourglass viscous (was refused in v1)."""
    _build_common(E=1000.0, nu=0.3)
    ops.element("LadrunoBrick", 1, *_CONN, 1, "-formulation", "uri",
                "-hourglass", "viscous", "-lumped")
    tags = ops.getEleTags() or []
    if isinstance(tags, int):
        tags = [tags]
    assert 1 in tags, "uri -hourglass viscous should be accepted"


def _viscous_bar_max_tip(hg_coeff, nsteps=400, dt=1.0e-4):
    """Explicit run of a clamped 4x1x1 uri bar under a sustained transverse tip
    load, -hourglass viscous with the given coefficient. Returns the peak |tip-z|
    over the run. In viscous mode there is NO hourglass stiffness, so the
    transverse (bending) motion of a one-element-deep bar is a damped zero-
    stiffness mode: it ramps toward a terminal velocity set by the viscous
    coefficient, so a LARGER coefficient ⇒ MORE damping ⇒ SMALLER peak tip-z. A
    wrong sign (anti-damping) would instead blow up. This makes the coefficient
    observable in the response — unlike a bare finiteness check."""
    nx, L = 4, 4.0
    E, nu, rho = 1000.0, 0.3, 1.0
    dx = L / nx

    def nt(i, j, k):
        return 1 + i + (nx + 1) * j + (nx + 1) * 2 * k

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k in range(2):
        for j in range(2):
            for i in range(nx + 1):
                ops.node(nt(i, j, k), i * dx, j * 1.0, k * 1.0)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)

    e = 1
    for i in range(nx):
        ops.element("LadrunoBrick", e,
                    nt(i, 0, 0), nt(i + 1, 0, 0), nt(i + 1, 1, 0), nt(i, 1, 0),
                    nt(i, 0, 1), nt(i + 1, 0, 1), nt(i + 1, 1, 1), nt(i, 1, 1),
                    1, "-formulation", "uri", "-hourglass", "viscous", hg_coeff,
                    "-lumped")
        e += 1

    for k in range(2):                          # clamp the x=0 face
        for j in range(2):
            ops.fix(nt(0, j, k), 1, 1, 1)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    tip = [nt(nx, j, k) for k in range(2) for j in range(2)]
    for n in tip:
        ops.load(n, 0.0, 0.0, 1.0 / len(tip))   # transverse tip load (excites bending=hg)

    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno")
    ops.analysis("Transient")

    maxd = 0.0
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        d = ops.nodeDisp(tip[0])[2]
        assert math.isfinite(d), "viscous explicit run produced a non-finite disp"
        maxd = max(maxd, abs(d))
    return maxd


def test_uri_viscous_hourglass_damps_and_is_stable():
    """Differential gate (NOT just does-not-crash): a larger viscous coefficient
    must produce a strictly smaller peak hourglass response, and both runs stay
    bounded. This fails if the coefficient never reaches the kernel (parser
    regression), if the damping sign is wrong (would blow up), or if the viscous
    branch is inert."""
    low = _viscous_bar_max_tip(0.05)
    high = _viscous_bar_max_tip(2.0)
    assert math.isfinite(low) and math.isfinite(high)
    assert low < 1.0e3 and high < 1.0e3, f"explicit run blew up (low={low:.3e}, high={high:.3e})"
    assert high < low, (
        f"more viscous damping should reduce the peak hourglass tip-z, but "
        f"coeff=2.0 gave {high:.3e} >= coeff=0.05 gave {low:.3e} "
        f"(coefficient not reaching the kernel, or wrong damping sign)"
    )


def test_hourglass_coefficient_reaches_kernel():
    """Regression for the factory coeff parse: a NUMERIC -hourglass coefficient
    (passed as a Python number, the idiomatic openseespy call) must actually
    reach the kernel. A larger FB `stiffness` hourglass coefficient stiffens the
    hourglass modes, so the static response under the mixed loads must be
    strictly smaller. If the coeff were dropped (e.g. OPS_GetString/strtod can't
    read a Python float), both runs would use the same default and tie."""
    soft = _solve_hex("LadrunoBrick", ["-formulation", "uri", "-hourglass", "stiffness", 0.02])
    stiff = _solve_hex("LadrunoBrick", ["-formulation", "uri", "-hourglass", "stiffness", 5.0])
    soft_mag = math.sqrt(sum(u * u for u in soft[0]))
    stiff_mag = math.sqrt(sum(u * u for u in stiff[0]))
    assert stiff_mag < soft_mag, (
        f"a larger hourglass coeff should stiffen the response, but "
        f"coeff=5.0 |u|={stiff_mag:.4e} >= coeff=0.02 |u|={soft_mag:.4e} "
        f"(numeric -hourglass coefficient not reaching the kernel)"
    )


@pytest.mark.parametrize("nu", [0.3, 0.499])
def test_eas_matches_sspbrick_distorted_hex(nu):
    """eas vs SSPbrick on the deliberately DISTORTED hex, component-by-component
    (disp + element force), at compressible and near-incompressible nu. Exercises
    the J[4..19] distortion terms of the SSPbrick port that the axis-aligned
    cantilever oracle (test_ladrunoBrick_bending.py) cannot reach — a regular
    grid has a diagonal Jacobian, zeroing every off-diagonal cross-block. eas
    condenses from the initial tangent, so for an elastic material the assembled
    stiffness equals SSPbrick's and the solved state must agree to ~1e-6."""
    ref = _solve_hex("SSPbrick", [], nu=nu)
    eas = _solve_hex("LadrunoBrick", ["-formulation", "eas"], nu=nu)
    _assert_close(eas[0], ref[0], f"disp[nu={nu}]", rtol=1e-6, atol=1e-8)
    _assert_close(eas[2], ref[2], f"force[nu={nu}]", rtol=1e-6, atol=1e-8)
