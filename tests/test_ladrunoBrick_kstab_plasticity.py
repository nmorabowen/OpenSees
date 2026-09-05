"""LadrunoBrick — Tier-A `Kstab` scope under a NON-softening material.

Pins the outcome of the 2026-07-30 study recorded in
`Ladruno_implementation/11_brick_asdconcrete_integration.md` §3.1, which asked
whether the damage-scaled stabilization should be generalized to a tangent-based
rule so that plasticity models are "covered" too.

The four facts pinned here:

  1. `damageScale() == 1` for a material with no "damage" channel is DELIBERATE:
     with upstream `J2Plasticity`, every formulation reproduces the exact
     homogeneous uniaxial solution and `ssp` stores ZERO stabilization energy.
     Nothing to degrade, so degrading would be wrong.
  2. The claimed mechanism is real: in a plastifying bending state the `ssp`
     stabilization's share of the internal energy GROWS (it does not stay put).
  3. ...but it converges away under mesh refinement (~O(h)), so it is a
     discretization quantity, not a persistent bias — no code change warranted.
  4. Tangent-tracking does NOT fix it: `uri`+`stiffness` already sizes its
     hourglass modulus from the CURRENT tangent shear `dd(3,3)`, and it records a
     LARGER stabilization share than the frozen-`Kstab` `ssp` at every mesh.

Plus the sharper trap the study surfaced: with ONE element through the bending
depth a single-point formulation never plastifies at all.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, NU, SIG0 = 200000.0, 0.3, 250.0
K3, G3 = E / (3 * (1 - 2 * NU)), E / (2 * (1 + NU))
HMOD = 0.002 * E                       # near-perfect; H==0 makes J2's local solve hunt

FORMS = {
    "std":  ["-formulation", "std"],
    "bbar": ["-formulation", "bbar"],
    "ssp":  ["-formulation", "ssp"],
    "uri":  ["-formulation", "uri", "-hourglass", "stiffness"],
    "eas":  ["-formulation", "eas"],
}


def _grid(nx, ny, nz, dx, dy, dz):
    nid, c = {}, 0
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                c += 1
                nid[(i, j, k)] = c
                ops.node(c, i * dx, j * dy, k * dz)
    return nid


def _conn(nid, i, j, k):
    return [nid[(i, j, k)], nid[(i + 1, j, k)], nid[(i + 1, j + 1, k)], nid[(i, j + 1, k)],
            nid[(i, j, k + 1)], nid[(i + 1, j, k + 1)], nid[(i + 1, j + 1, k + 1)],
            nid[(i, j + 1, k + 1)]]


def _solve(top, ctrl_dof, du, nsteps, eles):
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-9, 60, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("DisplacementControl", top[0], ctrl_dof, du)
    ops.analysis("Static")
    hist, wext, fprev, uprev = [], 0.0, 0.0, 0.0
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        P = len(top) * ops.getLoadFactor(1)
        u = ops.nodeDisp(top[0], ctrl_dof)
        wext += 0.5 * (P + fprev) * (u - uprev)
        uprev, fprev = u, P
        ehg = 0.0
        for e in eles:
            ops.eleResponse(e, "forces")
            ehg += ops.eleResponse(e, "hourglassEnergy")[0]
        hist.append((u, P, ehg, wext))
    return hist


# ---------------------------------------------------------------------------
# 1. Homogeneous plastic field: nothing for a degradation rule to do
# ---------------------------------------------------------------------------
LX, LY, LZ = 4.0, 1.0, 6.0


def _bar(form, nx=4, nz=6, utop=0.06, nsteps=30):
    """Roller-sided prismatic bar. Lateral contraction is free, so the exact
    solution is homogeneous uniaxial tension for ANY consistent element."""
    dx, dy, dz = LX / nx, LY, LZ / nz
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    nid = _grid(nx, 1, nz, dx, dy, dz)
    ops.nDMaterial("J2Plasticity", 1, K3, G3, SIG0, SIG0, 0.0, HMOD)
    eles, t = [], 0
    for k in range(nz):
        for i in range(nx):
            t += 1
            ops.element("LadrunoBrick", t, *_conn(nid, i, 0, k), 1, *FORMS[form])
            eles.append(t)
    for i in range(nx + 1):
        for j in (0, 1):
            ops.fix(nid[(i, j, 0)], 0, 0, 1)
    for j in (0, 1):
        for k in range(nz + 1):
            ops.fix(nid[(0, j, k)], 1, 0, 0)
    for i in range(nx + 1):
        for k in range(nz + 1):
            ops.fix(nid[(i, 0, k)], 0, 1, 0)
    top = [nid[(i, j, nz)] for i in range(nx + 1) for j in (0, 1)]
    for n in top[1:]:
        ops.equalDOF(top[0], n, 3)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    tot = 0.0
    for n in top:
        cx, cy = ops.nodeCoord(n)[0], ops.nodeCoord(n)[1]
        w = (0.5 if cx in (0.0, LX) else 1.0) * (0.5 if cy in (0.0, LY) else 1.0)
        ops.load(n, 0.0, 0.0, w)
        tot += w
    h = _solve(top, 3, utop / nsteps, nsteps, eles)
    # _solve normalises by len(top); rescale to the true tributary resultant
    return [(u, P / len(top) * tot, ehg, w) for (u, P, ehg, w) in h]


def _exact(u):
    eps = u / LZ
    if eps <= SIG0 / E:
        return E * eps * LX * LY
    epsp = (E * eps - SIG0) / (E + HMOD)
    return (SIG0 + HMOD * epsp) * LX * LY


@pytest.mark.parametrize("form", list(FORMS))
def test_homogeneous_plastic_field_is_exact_and_unstabilized(form):
    """A fully plastified HOMOGENEOUS field is reproduced exactly by every
    formulation, and `ssp` stores no stabilization energy at all. This is why
    damageScale()==1 for a non-softening material is correct rather than a
    missed case: there is no spurious stabilization force to degrade."""
    h = _bar(form)
    assert len(h) > 20, f"{form}: analysis stalled after {len(h)} steps"
    u, P, ehg, _ = h[-1]
    assert u / LZ > 2.0 * SIG0 / E, "did not reach a plastic state"
    assert P == pytest.approx(_exact(u), rel=2e-3)
    if form in ("ssp", "uri"):
        # hourglass content of an exactly homogeneous field is zero
        assert abs(ehg) < 1e-9 * abs(P) * u


# ---------------------------------------------------------------------------
# 2-4. Plastic bending: the mechanism is real, converges away, and is not
#      cured by tangent-tracking
# ---------------------------------------------------------------------------
DEPTH, THICK, HEIGHT = 1.0, 0.25, 6.0


def _cant(form, nd, elastic=False, utip=0.25, nsteps=100):
    nx, nz = nd, 3 * nd
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    nid = _grid(nx, 1, nz, DEPTH / nx, THICK, HEIGHT / nz)
    if elastic:
        ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    else:
        ops.nDMaterial("J2Plasticity", 1, K3, G3, SIG0, SIG0, 0.0, HMOD)
    eles, t = [], 0
    for k in range(nz):
        for i in range(nx):
            t += 1
            ops.element("LadrunoBrick", t, *_conn(nid, i, 0, k), 1, *FORMS[form])
            eles.append(t)
    for i in range(nx + 1):
        for j in (0, 1):
            ops.fix(nid[(i, j, 0)], 1, 1, 1)
    top = [nid[(i, j, nz)] for i in range(nx + 1) for j in (0, 1)]
    for n in top[1:]:
        ops.equalDOF(top[0], n, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in top:
        ops.load(n, 1.0, 0.0, 0.0)
    return _solve(top, 1, utip / nsteps, nsteps, eles)


def _hgfrac(h, frac):
    i = max(0, int(frac * len(h)) - 1)
    return h[i][2] / h[i][3]


def test_stabilization_share_grows_as_the_section_plastifies():
    """The challenged mechanism IS real: with a frozen elastic Kstab and a
    material that publishes no damage channel, the stabilization's share of the
    internal energy climbs as the membrane plastifies."""
    h = _cant("ssp", 4)
    assert len(h) > 50, f"analysis stalled after {len(h)} steps"
    early, late = _hgfrac(h, 0.15), _hgfrac(h, 1.0)
    assert early < 0.10, f"elastic-range share unexpectedly high: {early:.1%}"
    assert late > 1.5 * early, f"share did not grow: {early:.1%} -> {late:.1%}"


def test_stabilization_share_converges_away_under_refinement():
    """...but the magnitude is a DISCRETIZATION quantity: refining the bending
    depth drives both the energy share and the load bias vs full integration
    toward zero, which is why no code change was made.

    NB the load bias (`ssp/bbar`, elastic-calibrated) is the honest error measure.
    The energy share is a TREND diagnostic only — it is not a bound on the force
    error, which exceeds it ~2.5-3x at every mesh here.

    All meshes are compared at a COMMON step index, not at each run's own final
    step: every run uses the same `du` under DisplacementControl, so a common
    index is a common tip displacement. Comparing at per-run final indices would
    let an early-truncating fine mesh pass the monotonicity assert for the wrong
    reason, since the share grows with plastification (see the test above).
    """
    runs = {}
    for nd in (2, 4, 8):
        for form in ("ssp", "bbar"):
            for el in (False, True):
                r = _cant(form, nd, elastic=el)
                assert len(r) > 40, f"nd={nd} {form} elastic={el}: stalled at {len(r)}"
                runs[(nd, form, el)] = r
    i = min(len(r) for r in runs.values()) - 1     # common step -> common u_tip
    shares, bias = {}, {}
    for nd in (2, 4, 8):
        hp = runs[(nd, "ssp", False)]
        shares[nd] = hp[i][2] / hp[i][3]
        # divide out the purely elastic discretisation difference
        bias[nd] = ((hp[i][1] / runs[(nd, "bbar", False)][i][1])
                    / (runs[(nd, "ssp", True)][i][1] / runs[(nd, "bbar", True)][i][1]))
    assert shares[2] > shares[4] > shares[8], shares
    assert shares[8] < 0.5 * shares[4], shares
    assert bias[2] > bias[4] > bias[8], bias
    assert bias[8] < 1.25, bias
    # the coarse-mesh bias is LARGE and must stay visible: it drives the meshing
    # rule in ADR section 3.1 (production RC walls are 2-3 elements thick).
    assert bias[2] > 1.4, bias


def test_tangent_tracking_does_not_reduce_the_stabilization_share():
    """`uri`+`stiffness` sizes its hourglass modulus from the CURRENT tangent
    shear dd(3,3), and still records a LARGER share than frozen-Kstab `ssp` at
    every mesh.

    ILLUSTRATIVE, NOT PROBATIVE — and pinned as such so nobody cites it as proof.
    The comparison is confounded twice: `uri` tracks the single shear entry (a
    rule ADR section 3.1 actually endorses) rather than the challenged Frobenius
    norm, and its Flanagan-Belytschko kappa = 0.05*G*vol*sum(b^2) is a
    different-magnitude operator than `ssp`'s condensed Kstab, so much of the gap
    measures operator size rather than the failure of tracking. The decisive
    argument is analytic (section 3.1): in a bending/axial-flow state neither
    dd(3,3) nor ||D|| degrades meaningfully, so no tangent rule moves the answer.
    """
    for nd in (2, 4, 8):
        hs, hu = _cant("ssp", nd), _cant("uri", nd)
        i = min(len(hs), len(hu)) - 1          # common step -> common u_tip
        assert i > 40, f"nd={nd}: stalled at {i + 1} steps"
        s, u = hs[i][2] / hs[i][3], hu[i][2] / hu[i][3]
        assert u > s, f"nd={nd}: uri {u:.1%} should exceed ssp {s:.1%}"


def test_eas_recondensation_removes_the_plastic_bias_that_ssp_carries():
    """The clean frozen-vs-recondensed A/B, using only shipped code.

    `eas` (true Simo-Rifai) rebuilds Kaa/Kda/Kad from `getTangent()` — the CURRENT
    tangent — at all 8 live GPs and re-condenses K* = Kdd - Kda Kaa^-1 Kad on every
    assembly (`formEAStrue`). `ssp` condenses ONCE from C(0) and freezes it
    (`buildSSP`). The two are elastically near-identical on this mesh family, so
    their plastic divergence isolates the stabilization treatment itself.

    This is the experiment that decides the design question: a tangent-consistent
    stabilization IS worth a lot (ssp's +45..+77% bias collapses to ~+1..+8%), but
    the way to get it is full RE-CONDENSATION, not a scalar scale on a frozen
    operator — which is why no `-hgScale`-style option was added to `ssp`.
    """
    for nd in (2, 4):
        r = {(f, el): _cant(f, nd, elastic=el)
             for f in ("ssp", "eas", "bbar") for el in (False, True)}
        i = min(len(v) for v in r.values()) - 1     # common step -> common u_tip
        assert i > 40, f"nd={nd}: stalled at {i + 1} steps"

        def bias(f):
            return ((r[(f, False)][i][1] / r[("bbar", False)][i][1])
                    / (r[(f, True)][i][1] / r[("bbar", True)][i][1]))

        # elastic responses agree, so the plastic gap is not a discretisation gap
        el_ssp = r[("ssp", True)][i][1] / r[("bbar", True)][i][1]
        el_eas = r[("eas", True)][i][1] / r[("bbar", True)][i][1]
        assert el_eas == pytest.approx(el_ssp, rel=5e-3), (el_eas, el_ssp)

        b_ssp, b_eas = bias("ssp"), bias("eas")
        assert b_ssp > 1.30, f"nd={nd}: expected a large ssp bias, got {b_ssp:.3f}"
        assert b_eas < 1.15, f"nd={nd}: expected eas near-unbiased, got {b_eas:.3f}"
        assert b_eas < 0.5 * (b_ssp - 1.0) + 1.0, (b_ssp, b_eas)


def test_ssp_one_element_deep_never_plastifies_in_bending():
    """The sharper trap: with ONE element through the bending depth, `ssp`'s only
    material evaluation point — the mean-dilatation core at the centroid — sees
    ~zero bending strain, so it returns an exactly ELASTIC pushover under a
    plastic material while `std`/`bbar` form a hinge on the same mesh. No
    Kstab-degradation rule of any kind can fix this; a damage model would not
    damage at that centroid either.

    Scoped to `ssp` deliberately. `uri`+`stiffness` strains its centroid from the
    plain centroid B rather than the mean-dilatation core, so it does pick up
    transverse shear and DOES yield — at an aspect-ratio-dependent point — which
    is why it is asserted separately below.
    """
    pl = _cant("ssp", 1, elastic=False, utip=0.30, nsteps=40)
    el = _cant("ssp", 1, elastic=True, utip=0.30, nsteps=40)
    n = min(len(pl), len(el))
    assert n > 20
    assert pl[n - 1][1] == pytest.approx(el[n - 1][1], rel=1e-6), \
        "ssp: expected an exactly elastic response one element deep"
    # essentially the whole response is stabilization energy
    assert _hgfrac(pl, 1.0) > 0.80
    # the full-integration forms DO plastify on the same mesh
    for form in ("std", "bbar"):
        p2 = _cant(form, 1, elastic=False, utip=0.30, nsteps=40)
        e2 = _cant(form, 1, elastic=True, utip=0.30, nsteps=40)
        m = min(len(p2), len(e2))
        assert p2[m - 1][1] < 0.5 * e2[m - 1][1], f"{form}: expected plastification"


def test_uri_one_element_deep_does_yield():
    """Companion to the above — pins that the trap is NOT generic to every
    single-point formulation, so nobody widens the guidance on a false premise."""
    pl = _cant("uri", 1, elastic=False, utip=0.30, nsteps=40)
    el = _cant("uri", 1, elastic=True, utip=0.30, nsteps=40)
    n = min(len(pl), len(el))
    assert n > 20
    assert pl[n - 1][1] < 0.5 * el[n - 1][1], "uri: expected plastification"
