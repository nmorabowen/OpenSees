"""ADR-41 C3.1 — mortar friction FORCE (Coulomb/Tresca, penalty λ_T≡0): the load-bearing phase.

The `MORTAR` LadrunoContactFE mode now adds the tangential traction. For each in-contact slave
node it builds the LOCAL weighted tangential gap, runs the SHIPPED LadrunoFrictionKernel return
map (unified cone cap = min(μN+c, τmax)) with the nodal normal pressure N_I = −p_normal,I (the C2
value), and scatters the (already-negated, motion-opposing) traction via the D/M operators:
F^s_K += Σ_I D_KI tFric_I, F^m_L += −Σ_I M_IL tFric_I. The per-global-slave-node elastic slip
gpT + engagement origin gT0 live on LadrunoContactDomain::MortarNormalState (the C3 analogue of
λ_N; committed in Domain::commit). λ_T (tangential Uzawa) is C3.3; the friction TANGENT is C3.2,
so these gates are EXPLICIT (CentralDifferenceLadruno force-only, like the NTS P3 ship).

Gates (capstone C3 / ADR-41):
  - a driven block slips at a = (Q − cap·A)/m — Tresca cap (μ=0) so the magnitude is independent
    of the normal pressure (the C2.2 Uzawa augments λ_N each commit; a Coulomb cap would drift,
    a Tresca cap does not), through the REAL CDL + D/M assembly;
  - the friction OPPOSES the motion (a Coulomb-loaded slave given a tangential velocity DECELERATES,
    never accelerates — the A1 BLOCKER-1 sign gate through real assembly);
  - μ=0 ∧ c=0 ∧ τmax=0 ⇒ BYTE-IDENTICAL to the frictionless C2 mortar path (the short-circuit).

The mechanics were validated oracle-first in contact_prototypes/proto_c3_mortar_friction.py (20/20).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

EPS = 1.0e4      # mortar normal penalty
DELTA = 1.0e-3   # seeded penetration (keeps the facet engaged ⇒ friction active)


def _master_quad(z=0.0, tags=(1, 2, 3, 4)):
    c = [(tags[0], 0.0, 0.0, z), (tags[1], 1.0, 0.0, z),
         (tags[2], 1.0, 1.0, z), (tags[3], 0.0, 1.0, z)]
    for t, x, y, zz in c:
        ops.node(t, x, y, zz)
        ops.fix(t, 1, 1, 1)
    return list(tags)


def _slave_quad(z, tags=(5, 6, 7, 8), m_node=0.25, fix_y=True, fix_z=True):
    c = [(tags[0], 0.0, 0.0, z), (tags[1], 1.0, 0.0, z),
         (tags[2], 1.0, 1.0, z), (tags[3], 0.0, 1.0, z)]
    for t, x, y, zz in c:
        ops.node(t, x, y, zz)
        ops.mass(t, m_node, m_node, m_node)
        # x free; y/z fixed so the only dynamics is tangential x (z fixed ⇒ the normal penetration
        # stays = DELTA ⇒ the node stays in contact ⇒ friction active, decoupled from x).
        ops.fix(t, 0, 1 if fix_y else 0, 1 if fix_z else 0)
    return list(tags)


def _run_cdl(nsteps, dt, sample):
    rec = []
    for k in range(nsteps):
        assert ops.analyze(1, dt) == 0, f"CDL step {k} failed"
        rec.append(sample(k))
    return rec


def _driven_mortar(mu, cohesion, Q, tauMax=0.0, nsteps=160, dt=1.0e-3, k1=50, k2=150):
    """A slave quad facet (z fixed at −DELTA ⇒ in contact) driven by a constant tangential +x
    force Q (per facet, split over the 4 nodes), on a fixed master quad. Returns the measured
    x-acceleration of the slave (slope of vx over steps [k1,k2])."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master_quad(0.0)
    _slave_quad(-DELTA)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    opts = ["-mortar", "-epsN", EPS, "-outward", 0.0, 0.0, 1.0]
    if mu > 0.0:
        opts += ["-mu", mu, "-epsT", EPS]
    if cohesion > 0.0:
        opts += ["-cohesion", cohesion, "-epsT", EPS]
    if tauMax > 0.0:
        opts += ["-tauMax", tauMax]
    ops.contact(1, 1, 2, *opts)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, Q / 4.0, 0.0, 0.0)             # tangential drive (+x), split over the facet
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    rec = _run_cdl(k2 + 1, dt, lambda k: ops.nodeVel(5, 1))
    return (rec[k2] - rec[k1]) / ((k2 - k1) * dt)


def test_c3_1_tresca_driven_block_slips_at_capped_accel():
    """A driven slave facet slips at a = (Q − cap·A)/m with a TRESCA cap (cap = cohesion,
    pressure-INDEPENDENT) — so the C2.2 normal-Uzawa growth of λ_N does NOT perturb the friction
    magnitude. Frictionless (cap=0) gives a = Q/m. The gap between them IS the friction force,
    assembled through the real mortar D/M operators + CDL."""
    Q, m_total, A = 10.0, 1.0, 1.0
    cohesion = 2.0                                  # Tresca cap = 2.0 (cap·A = 2.0)
    a_fric = _driven_mortar(mu=0.0, cohesion=cohesion, Q=Q)
    a_free = _driven_mortar(mu=0.0, cohesion=0.0, Q=Q)
    a_th_fric = (Q - cohesion * A) / m_total        # 8.0
    a_th_free = Q / m_total                          # 10.0
    assert abs(a_free - a_th_free) / a_th_free < 0.03, f"frictionless a={a_free} vs {a_th_free}"
    assert abs(a_fric - a_th_fric) / a_th_fric < 0.03, f"Tresca a={a_fric} vs {a_th_fric}"
    # the friction decelerated the block by exactly cap·A/m (the assembled tangential force).
    assert a_free - a_fric == pytest.approx(cohesion * A / m_total, rel=0.05)


def test_c3_1_friction_opposes_motion_sign_gate():
    """THE SIGN GATE (A1 BLOCKER-1, through real CDL + D/M assembly): a Coulomb-loaded slave given
    a tangential +x velocity must DECELERATE (friction opposes motion) — never accelerate (a wrong
    sign injects energy). N grows under the normal Uzawa, but the DIRECTION (opposing) is invariant."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master_quad(0.0)
    _slave_quad(-DELTA)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-mu", 0.3, "-epsT", EPS,
                "-outward", 0.0, 0.0, 1.0)
    for t in (5, 6, 7, 8):
        ops.setNodeVel(t, 1, 1.0, "-commit")        # initial +x velocity
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    v_prev = ops.nodeVel(5, 1)
    assert v_prev > 0.0
    for k in range(30):
        assert ops.analyze(1, 1.0e-4) == 0
        v = ops.nodeVel(5, 1)
        # monotone non-increasing while still moving +x (friction OPPOSES, never speeds it up).
        assert v <= v_prev + 1.0e-9, f"friction ACCELERATED the slave (sign bug): {v_prev}->{v}"
        v_prev = v
    assert v_prev < 1.0, "slave did not decelerate at all (friction inert?)"


def test_c3_1_frictionless_byte_identical():
    """μ=0 ∧ c=0 ∧ τmax=0 ⇒ the friction block short-circuits before any state touch, so the model
    evolves BITWISE identically to a mortar contact defined with NO friction flags (the C2 path)."""
    def _run(with_zero_friction):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        _master_quad(0.0)
        _slave_quad(-DELTA, fix_z=False)            # z free ⇒ real normal dynamics to compare
        ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
        ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
        if with_zero_friction:
            ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-mu", 0.0, "-cohesion", 0.0,
                        "-outward", 0.0, 0.0, 1.0)
        else:
            ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-outward", 0.0, 0.0, 1.0)
        for t in (5, 6, 7, 8):
            ops.setNodeVel(t, 1, 0.3, "-commit")    # some tangential motion (would trip friction)
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.integrator("CentralDifferenceLadruno")
        ops.algorithm("Linear")
        ops.analysis("Transient")
        out = []
        for _ in range(12):
            assert ops.analyze(1, 1.0e-4) == 0
            out.append(tuple(ops.nodeDisp(t, d) for t in (5, 6, 7, 8) for d in (1, 3)))
        return out
    assert _run(True) == _run(False), "zero-friction mortar is not byte-identical to frictionless"
