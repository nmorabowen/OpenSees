"""LadrunoBrick — cantilever bending-convergence benchmark (Belytschko Table 8.2
style). This is the oracle that *validates an assumed-strain element*: the patch
and rank tests cannot (gamma-orthogonality makes any hourglass variant pass
them), but bending convergence catches a wrong stabilization operator.

Setup: a slender cantilever (L=10, 1x1 square section) meshed with nx hexes
along the length, ONE element through the depth — the harsh shear-locking case
for trilinear hexes. nu=0 isolates shear locking (no volumetric component).
Transverse tip load P=1; analytic Euler tip deflection delta = P L^3/(3 E I),
I = b h^3/12 = 1/12 -> delta = 4.0.

Calibrated behaviour (2026-06-01, this build):
    mesh    std     bbar    uri(stiffness)
    nx=4    0.244   0.270   0.127
    nx=8    0.564   0.736   0.407
    nx=16   0.841   1.296   0.940
std shear-locks; uri+stiffness is *worse* at coarse mesh (bending lives in the
hourglass space, governed by the perturbation control). The assumed-strain
`physical` formulation should reach ~0.9 even at nx=4 — the validating contrast.

Plan ref: Ladruno_implementation/09_ladruno_brick.md; theory Belytschko sec 8.7.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

ANALYTIC = 4.0   # Euler-Bernoulli tip deflection, P L^3 / (3 E I), I=1/12


def _cantilever_tip(form_args, nx, ny=1, nz=1, L=10.0, b=1.0, h=1.0,
                    E=1000.0, nu=0.0, P=1.0, elem="LadrunoBrick"):
    """Tip deflection (z) of a clamped cantilever meshed nx x ny x nz, tip
    transverse (z) load P, built from `elem` (LadrunoBrick of the given
    formulation, or e.g. SSPbrick for a reference assumed-strain hex)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    dx, dy, dz = L / nx, b / ny, h / nz

    def nt(i, j, k):
        return 1 + i + (nx + 1) * j + (nx + 1) * (ny + 1) * k

    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(nt(i, j, k), i * dx, j * dy, k * dz)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)

    e = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ops.element(elem, e,
                            nt(i, j, k), nt(i + 1, j, k), nt(i + 1, j + 1, k), nt(i, j + 1, k),
                            nt(i, j, k + 1), nt(i + 1, j, k + 1), nt(i + 1, j + 1, k + 1), nt(i, j + 1, k + 1),
                            1, *form_args)
                e += 1

    for k in range(nz + 1):                 # clamp the x=0 face
        for j in range(ny + 1):
            ops.fix(nt(0, j, k), 1, 1, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    tip = [nt(nx, j, k) for k in range(nz + 1) for j in range(ny + 1)]
    for n in tip:
        ops.load(n, 0.0, 0.0, P / len(tip))

    ops.system("FullGeneral")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return sum(ops.nodeDisp(n)[2] for n in tip) / len(tip)


def _ratio(form_args, nx):
    return _cantilever_tip(form_args, nx) / ANALYTIC


# --------------------------------------------------------------------------
# the locking baseline (documents the problem the assumed strain solves)
# --------------------------------------------------------------------------
def test_std_shear_locks_in_bending():
    """Full integration shear-locks badly with one element through the depth."""
    assert _ratio(["-formulation", "std"], 4) < 0.40
    # ... and converges slowly from below as the mesh refines
    assert _ratio(["-formulation", "std"], 4) < _ratio(["-formulation", "std"], 16) < 1.0


# --------------------------------------------------------------------------
# the validating gate for physical (assumed strain)
# --------------------------------------------------------------------------
def test_physical_relieves_shear_locking():
    """Belytschko-Bindeman assumed-strain stabilization: near-exact bending even
    on a coarse mesh, dramatically better than std/uri-stiffness."""
    phys4 = _ratio(["-formulation", "uri", "-hourglass", "physical"], 4)
    std4 = _ratio(["-formulation", "std"], 4)
    uri4 = _ratio(["-formulation", "uri", "-hourglass", "stiffness"], 4)
    # assumed strain should recover most of the bending even at nx=4
    assert phys4 > 0.80, f"physical nx=4 ratio {phys4:.3f} should be ~0.9"
    # and crush both the full-integration and perturbation-stabilized results
    assert phys4 > 2.0 * std4, f"physical {phys4:.3f} should beat std {std4:.3f}"
    assert phys4 > 3.0 * uri4, f"physical {phys4:.3f} should beat uri-stiffness {uri4:.3f}"


def test_physical_converges():
    """physical stays accurate (not overshooting) under refinement at nu=0."""
    for nx in (4, 8, 16):
        r = _ratio(["-formulation", "uri", "-hourglass", "physical"], nx)
        assert 0.80 < r < 1.15, f"physical nx={nx} ratio {r:.3f} out of [0.80,1.15]"


def test_physical_matches_assumed_strain_reference():
    """Cross-validate against OpenSees' proven assumed-strain hex (SSPbrick):
    in compressible bending the two should agree closely."""
    phys = _cantilever_tip(["-formulation", "uri", "-hourglass", "physical"], 8)
    ref = _cantilever_tip([], 8, elem="SSPbrick")
    assert abs(phys - ref) / ref < 0.05, (
        f"physical tip {phys:.4f} should match SSPbrick {ref:.4f} within 5%"
    )


def test_physical_volumetric_locking_is_a_known_limitation():
    """DOCUMENTED LIMITATION: physical cures SHEAR locking but NOT volumetric
    locking — near-incompressible it locks (use -formulation bbar there). This
    test pins that behaviour so a future volumetric fix is detected as a change.
    A correct general element (SSPbrick) stays accurate at nu->0.5."""
    nu = 0.499
    phys = _cantilever_tip(["-formulation", "uri", "-hourglass", "physical"], 8, nu=nu) / ANALYTIC
    ref = _cantilever_tip([], 8, nu=nu, elem="SSPbrick") / ANALYTIC
    assert phys < 0.5, f"physical is expected to volumetrically lock at nu={nu} (got {phys:.3f})"
    assert ref > 0.9, f"SSPbrick (general) should not lock at nu={nu} (got {ref:.3f})"


# --------------------------------------------------------------------------
# ssp (v2) — bbar + statically-condensed enhanced strain (SSPbrick port).
# The GENERAL-nu element: cures BOTH shear and volumetric locking. This is the
# validating gate the patch/rank tests cannot provide.
# --------------------------------------------------------------------------
def test_ssp_converges_in_bending():
    """ssp relieves shear locking like physical: near-exact bending even on a
    coarse mesh, accurate (not overshooting) under refinement at nu=0."""
    for nx in (4, 8, 16):
        r = _ratio(["-formulation", "ssp"], nx)
        assert 0.80 < r < 1.15, f"ssp nx={nx} ratio {r:.3f} out of [0.80,1.15]"


@pytest.mark.parametrize("nu", [0.0, 0.3, 0.45, 0.499])
def test_ssp_matches_sspbrick_across_all_nu(nu):
    """THE validating gate. LadrunoBrick ssp is a port of SSPbrick's bbar +
    statically-condensed enhanced-strain structure, with the enhanced modes
    condensed from the INITIAL tangent. For a linear-elastic material that makes
    the assembled stiffness *identical* to SSPbrick, so the tip deflection must
    agree to near machine precision — across the FULL nu range, including the
    near-incompressible regime where the fixed-assumed-strain `physical`
    formulation volumetrically locks."""
    ssp = _cantilever_tip(["-formulation", "ssp"], 8, nu=nu)
    ref = _cantilever_tip([], 8, nu=nu, elem="SSPbrick")
    assert abs(ssp - ref) <= 1e-6 * abs(ref) + 1e-9, (
        f"ssp tip {ssp:.8f} should match SSPbrick {ref:.8f} at nu={nu} "
        f"(|d|={abs(ssp-ref):.2e})"
    )


def test_ssp_cures_volumetric_locking_where_physical_fails():
    """The headline reason ssp exists: at nu->0.5 the `physical` shear-cure
    volumetrically locks (<0.5) but ssp stays accurate (>0.9), matching the
    proven general element SSPbrick. ssp wins on the exact axis physical loses."""
    nu = 0.499
    ssp = _cantilever_tip(["-formulation", "ssp"], 8, nu=nu) / ANALYTIC
    phys = _cantilever_tip(["-formulation", "uri", "-hourglass", "physical"], 8, nu=nu) / ANALYTIC
    assert ssp > 0.9, f"ssp should NOT volumetrically lock at nu={nu} (got {ssp:.3f})"
    assert phys < 0.5, f"physical is expected to lock at nu={nu} (got {phys:.3f})"
