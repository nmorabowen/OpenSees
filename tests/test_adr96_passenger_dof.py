"""ADR 96 -- contact and ZeroLength on ndf >= 3 nodes, the pressure DOF as a
passenger (`Ladruno_implementation/96_ladruno_contact_passenger_dof_adr.md`;
TIMs request 2026-09-07, F1).

Gates (the ADR's G2 / G3; G1 is the byte-identity harness + the shipped battery):

  G3  `zeroLength` + `ENTMaterial` between an ndf-3 node and an ndf-4 node opens
      under tension (zero spring force, the ndf-3 node moves freely) and carries
      compression (equilibrium: spring force == applied load), with the ndf-4
      node's 4th DOF untouched. The same deck on a (4,4) and a (3,3) pair gives
      the same mechanical answer; a (3,4) pair with a ROTATIONAL -dir is refused.

  G2  A `LadrunoUP` column (H8, bbar, equal-order p, drained static path) under an
      ndf-3 rigid platen in NTS contact on its top face, pushed statically:
      the summed contact normal traction equals the applied load, no negative
      contact pressure, and the pore-pressure field is identical to the same
      column bonded by `equalDOF` on DOFs 1-3.

Every deck is tiny and serial; zone_a.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

K_SPRING = 1.0e4
P_LOAD = 250.0
PASSENGER_VAL = 7.0     # imposed on DOF 4 of the ndf-4 node: never read, never written


# ---------------------------------------------------------------------------
#  G3 -- the no-tension spring across a (3,4) pair
# ---------------------------------------------------------------------------

def _spring_pair(ndf_a, ndf_b, direction=3, load=+P_LOAD, ent=True):
    """Node 1 (ndf_a) fixed; node 2 (ndf_b) coincident, free in `direction`,
    joined by a zeroLength spring in `direction`. Load node 2 along it.

    On an ndf > 3 node 2, DOF 4 is neither fixed nor loaded: it is IMPOSED
    (`sp`) at PASSENGER_VAL so that (i) the element demonstrably never READS
    it -- the spring answer is the same as on a (3,3) pair -- and (ii) never
    WRITES it -- the reaction the `sp` collects on DOF 4 is exactly zero.

    Returns (rc, disp of node 2 along `direction`, node-2 DOF-4 value or
    None, reaction on DOF 4 or None, spring force from the element response).
    """
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', ndf_a)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.fix(1, *([1] * ndf_a))
    ops.model('basic', '-ndm', 3, '-ndf', ndf_b)
    ops.node(2, 0.0, 0.0, 0.0)
    fx = [1] * ndf_b
    fx[direction - 1] = 0
    if ndf_b > 3:
        fx[3] = 0          # DOF 4 is imposed by `sp` below, not fixed
    ops.fix(2, *fx)
    if ent:
        ops.uniaxialMaterial('ENT', 1, K_SPRING)
    else:
        ops.uniaxialMaterial('Elastic', 1, K_SPRING)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', direction)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ld = [0.0] * ndf_b
    ld[direction - 1] = load
    ops.load(2, *ld)
    if ndf_b > 3:
        ops.sp(2, 4, PASSENGER_VAL)
    ops.constraints('Transformation')   # Plain silently homogenises the DOF-4 sp
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-12, 20, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.25)
    ops.analysis('Static')
    rc = 0
    for _ in range(4):
        rc = ops.analyze(1)
        if rc != 0:
            break
    u = ops.nodeDisp(2, direction)
    p4 = r4 = None
    if ndf_b > 3:
        p4 = ops.nodeDisp(2, 4)
        ops.reactions()
        r4 = ops.nodeReaction(2, 4)
    force = list(ops.eleResponse(1, 'force'))
    return rc, u, p4, r4, force


def test_g3_ent_spring_3_4_carries_compression():
    rc, u, p4, r4, force = _spring_pair(3, 4, load=-P_LOAD)
    assert rc == 0
    # compression: the spring closes, K*u balances the load
    assert abs(u + P_LOAD / K_SPRING) <= 1.0e-12 * P_LOAD / K_SPRING, u
    assert p4 == PASSENGER_VAL, p4               # the passenger sits where sp put it
    assert r4 == 0.0, r4                         # ... and the element never wrote to it
    # element-sized `force` response (3 + 4 slots): node 2's z slot carries
    # the load's sign and the passenger slot [6] is identically zero
    assert len(force) == 7, force
    assert abs(force[2] - P_LOAD) <= 1.0e-9 * P_LOAD, force
    assert abs(force[5] + P_LOAD) <= 1.0e-9 * P_LOAD, force
    assert force[6] == 0.0, force


def test_g3_ent_spring_3_4_opens_under_tension():
    """ENT under tension: zero stiffness ⇒ the singular free DOF is what a
    no-tension interface IS. The deck therefore fixes nothing extra and just
    asserts the spring transmits NOTHING: the element force is zero and the
    reaction at the fixed node is zero, on the first (elastic-predictor) step
    that ENT lets converge."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.model('basic', '-ndm', 3, '-ndf', 4)
    ops.node(2, 0.0, 0.0, 0.0)
    ops.fix(2, 1, 1, 0, 0)          # DOF 4 imposed by sp below
    # a soft parallel elastic spring keeps the tension DOF non-singular so the
    # opened ENT can be measured rather than inferred from a failed solve
    ops.uniaxialMaterial('ENT', 1, K_SPRING)
    ops.uniaxialMaterial('Elastic', 2, 1.0e-3 * K_SPRING)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 3)
    ops.element('zeroLength', 2, 1, 2, '-mat', 2, '-dir', 3)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 0.0, 0.0, +P_LOAD, 0.0)
    ops.sp(2, 4, PASSENGER_VAL)
    ops.constraints('Transformation')   # Plain silently homogenises the DOF-4 sp; ops.numberer('Plain'); ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-12, 30, 0); ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.25); ops.analysis('Static')
    for _ in range(4):
        assert ops.analyze(1) == 0
    u = ops.nodeDisp(2, 3)
    # only the soft spring carries: u = P / (1e-3 K)
    assert abs(u - P_LOAD / (1.0e-3 * K_SPRING)) <= 1.0e-9 * u, u
    ent = list(ops.eleResponse(1, 'force'))
    assert max(abs(v) for v in ent) == 0.0, ent      # ENT opened: nothing
    assert ops.nodeDisp(2, 4) == PASSENGER_VAL       # passenger where sp put it
    ops.reactions()
    assert ops.nodeReaction(2, 4) == 0.0             # ... and never written


@pytest.mark.parametrize('pair', [(3, 3), (4, 4), (4, 3), (3, 6)])
def test_g3_same_answer_on_every_pair(pair):
    rc, u, p4, r4, force = _spring_pair(*pair, load=-P_LOAD)
    assert rc == 0, pair
    assert abs(u + P_LOAD / K_SPRING) <= 1.0e-12 * P_LOAD / K_SPRING, (pair, u)
    if p4 is not None:
        assert p4 == PASSENGER_VAL and r4 == 0.0, (pair, p4, r4)
    # (3,3) is the vanilla D3N6 path, (3,6)/(4,3)/(4,4) the passenger path:
    # the `force` response is element-sized on every pair, node 2's z slot
    # sits at ndf_a + 2, and every non-translational slot is exactly zero
    na, nb = pair
    assert len(force) == na + nb, (pair, force)
    assert abs(force[2] - P_LOAD) <= 1.0e-9 * P_LOAD, (pair, force)
    assert abs(force[na + 2] + P_LOAD) <= 1.0e-9 * P_LOAD, (pair, force)
    extras = [force[i] for i in range(na + nb) if i not in (0, 1, 2, na, na + 1, na + 2)]
    assert all(v == 0.0 for v in extras), (pair, force)


def test_g3_rotational_dir_is_refused_on_a_passenger_pair(capfd):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    ops.node(1, 0.0, 0.0, 0.0); ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.model('basic', '-ndm', 3, '-ndf', 4)
    ops.node(2, 0.0, 0.0, 0.0); ops.fix(2, 1, 1, 1, 0)
    ops.uniaxialMaterial('Elastic', 1, K_SPRING)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 4)
    out = capfd.readouterr()
    text = out.err + out.out
    assert 'passenger mode, ADR-96' in text, text
    assert 'element disabled' in text


def test_g3_ndf_2_pair_is_still_refused_as_vanilla(capfd):
    """The relaxation is 3-D and ndf >= 3 only: a (2,3) pair still gets the
    vanilla 'differing dof at ends' refusal."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 2)
    ops.node(1, 0.0, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.node(2, 0.0, 0.0, 0.0); ops.fix(2, 1, 1, 0)
    ops.uniaxialMaterial('Elastic', 1, K_SPRING)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 3)
    text = ''.join(capfd.readouterr())
    assert 'differing dof at ends' in text, text


# ---------------------------------------------------------------------------
#  G2 -- NTS contact between an ndf-3 rigid platen and a LadrunoUP column
# ---------------------------------------------------------------------------
#
#  Deck (the ADR's G2 gate, the TIMs PM-01 shape):
#
#      z=4 . . . . 4 ndf-3 PLATEN nodes, coincident with the column top,
#                  lateral DOFs held, uz free, each loaded -P/4
#             ^^^  NTS penalty contact, kn = 1e3 * E*A/L
#      z=4  +----+ 4 ndf-4 column top nodes = ONE quad master segment
#           |    |
#           | UP | 4 x LadrunoUP H8 (bbar, equal-order p, ElasticIsotropic)
#           |    | ux = uy = 0 on every node (1-D confinement rollers)
#      z=0  +----+ uz fixed
#
#  The passenger under test is the column node's 4th DOF (p).  On the STATIC
#  drained path the u-p tangent is [K, -Q ; 0, H]: the p-row has no u
#  dependence, so with a Dirichlet p anywhere the pressure field is identically
#  zero and the column is exactly drained -- K u = P.  That is the reference
#  the twin (equalDOF on DOFs 1-3, no contact) must reproduce EXACTLY: the
#  contact adapter must never have read, written or coupled DOF 4.
#
#  Both drainage placements are run: `top` is the ADR's literal wording (p
#  FIXED on the four interface nodes) and `base` is the sharper one (p FREE on
#  the interface nodes, so a live pore-pressure equation number sits in the
#  DOF_Group slot that the pre-ADR-96 FE_Element::setID() copied into a
#  translation slot before returning -3).

UP_E, UP_NU, UP_RHO = 1.0e4, 0.25, 2.0        # kPa, Mg/m3
UP_KF, UP_PORO, UP_RHOF, UP_PERM = 2.2e5, 0.4, 1.0, 1.0e-4
COL_B, COL_L, COL_NZ = 1.0, 4.0, 4            # 1 x 1 x 4 of unit bricks
G2_P = 100.0                                  # total platen load (kN)
# kn = 1e3 * (E A / L) = 1e3 * 1e4 * 1 / 4 = 2.5e6  ->  per-slave penetration
# P/(4 kn) = 1.0e-5 m against a drained column shortening of
# P L / (Eoed A) = 100*4/1.2e4 = 3.3333e-2 m.  MEASURED at this kn (both
# drainage placements, build f001cd067):
#     sum(Fn)      = 100.000000000655       (rel err 6.6e-12 of P)
#     per-slave Fn = 25.0000000002 x4       (min > 0: compression only)
#     penetration  = 9.999000e-06 = P/(4 kn) - gap0   (exact to 1e-12 rel)
#     max|dp|      = 2.6e-24  (p is identically 0 on the drained static path)
#     max|du_z|    = 6.9e-18, i.e. 2.1e-16 RELATIVE -- the penalty acts only on
#                    the platen's own extra travel, never on the column, so the
#                    1e-6 gate below has ~10 orders of margin at this kn.
G2_KN = 1.0e3 * UP_E * COL_B * COL_B / COL_L
G2_GAP0 = 1.0e-9                              # start just penetrated (the
                                              # shipped ADR-39 NTS idiom)
PLATEN_BASE = 9001


def _g2_hex(nid, iz):
    """8 node tags of the iz-th hex, bottom face CCW then top face CCW."""
    return [nid[(0, 0, iz)], nid[(1, 0, iz)], nid[(1, 1, iz)], nid[(0, 1, iz)],
            nid[(0, 0, iz + 1)], nid[(1, 0, iz + 1)],
            nid[(1, 1, iz + 1)], nid[(0, 1, iz + 1)]]


def _g2_column(drain):
    """1x1x4 LadrunoUP H8 column on an ndf-4 model (caller opened it).

    ux = uy = 0 on every node (1-D confinement rollers), uz fixed at the base,
    p = 0 fixed on the drained face (`top` or `base`), p free elsewhere.
    Returns ({(ix,iy,iz): tag}, [4 top tags in CCW winding])."""
    nid = {}
    k = 1
    for iz in range(COL_NZ + 1):
        for iy in range(2):
            for ix in range(2):
                ops.node(k, ix * COL_B, iy * COL_B, iz * COL_L / COL_NZ)
                nid[(ix, iy, iz)] = k
                k += 1
    ops.nDMaterial('ElasticIsotropic', 1, UP_E, UP_NU, UP_RHO)
    for iz in range(COL_NZ):
        ops.element('LadrunoUP', iz + 1, *_g2_hex(nid, iz), 1,
                    '-Kf', UP_KF, '-poro', UP_PORO, '-rhoF', UP_RHOF,
                    '-perm', UP_PERM, UP_PERM, UP_PERM,
                    '-formulation', 'bbar', '-pOrder', 'equal',
                    '-stab', 'auto', '-dynSeepage', 'off')
    for (ix, iy, iz), n in nid.items():
        uz = 1 if iz == 0 else 0
        fp = 1 if ((drain == 'top' and iz == COL_NZ)
                   or (drain == 'base' and iz == 0)) else 0
        ops.fix(n, 1, 1, uz, fp)
    top = [nid[(0, 0, COL_NZ)], nid[(1, 0, COL_NZ)],
           nid[(1, 1, COL_NZ)], nid[(0, 1, COL_NZ)]]
    return nid, top


def _g2_platen(z, fix_lateral):
    """4 ndf-3 nodes over the column's top corners (the ADR's rigid platen).
    Opens an ndf-3 model on the SAME domain -- the mixed-ndf idiom of
    `_spring_pair` above."""
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    tags = []
    for i, (ix, iy) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)]):
        t = PLATEN_BASE + i
        ops.node(t, ix * COL_B, iy * COL_B, z)
        if fix_lateral:
            ops.fix(t, 1, 1, 0)      # rigid in-plane, uz free (frictionless
                                     # contact gives no lateral stiffness)
        tags.append(t)
    return tags


def _g2_run_contact(drain, kn=G2_KN, P=G2_P, gap0=G2_GAP0):
    """The contact model.  Returns (rc, nid, top, platen, per-slave normal
    force list)."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 4)
    nid, top = _g2_column(drain)
    platen = _g2_platen(COL_L - gap0, fix_lateral=True)
    ops.contactSurface(10, '-master', 4, *top)      # ndf-4 master facet
    ops.contactSurface(20, '-slave', *platen)       # ndf-3 slaves
    ops.contact(1, 10, 20, kn, 0.0, 0.0, '-outward', 0.0, 0.0, 1.0)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for t in platen:
        ops.load(t, 0.0, 0.0, -P / 4.0)
    ops.constraints('LadrunoContact')
    ops.numberer('Plain')
    ops.system('FullGeneral')          # the honest-p tangent is unsymmetric
    ops.integrator('LoadControl', 1.0)
    ops.test('NormDispIncr', 1.0e-12, 50, 0)
    ops.algorithm('Newton')
    ops.analysis('Static')
    rc = ops.analyze(1)
    fn = [ops.ladrunoContactForce(t) for t in platen]
    return rc, nid, top, platen, fn


def _g2_run_bonded(drain, P=G2_P):
    """The TWIN: the same column, the platen bonded to the top nodes by
    `equalDOF` on DOFs 1-3 (the retained node is the ndf-4 column node, so the
    platen inherits the column's roller SPs and no redundant SP is added --
    the rank-deficient-KKT trap)."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 4)
    nid, top = _g2_column(drain)
    platen = _g2_platen(COL_L, fix_lateral=False)
    for tn, pn in zip(top, platen):
        ops.equalDOF(tn, pn, 1, 2, 3)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for t in platen:
        ops.load(t, 0.0, 0.0, -P / 4.0)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.integrator('LoadControl', 1.0)
    ops.test('NormDispIncr', 1.0e-12, 50, 0)
    ops.algorithm('Newton')
    ops.analysis('Static')
    rc = ops.analyze(1)
    return rc, nid, top, platen


@pytest.mark.parametrize('drain', ['top', 'base'])
def test_g2_up_column_under_an_ndf3_platen(drain):
    """The ADR's G2 gate.  `drain='top'` is its literal wording (p fixed on the
    interface nodes); `drain='base'` leaves a LIVE pore-pressure equation on
    every interface node -- the exact DOF_Group slot the pre-ADR-96
    `FE_Element::setID()` copied into a translation slot before returning -3."""
    rc, nid, top, platen, fn = _g2_run_contact(drain)
    assert rc == 0, f'contact model did not converge (drain={drain})'

    # (2) equilibrium: the summed normal traction IS the applied load
    tot = sum(fn)
    assert abs(tot - G2_P) <= 1.0e-8 * G2_P, (drain, tot, fn)

    # (3) compression only -- no negative contact pressure
    assert min(fn) >= 0.0, (drain, fn)
    assert min(fn) > 0.0, f'a slave carries nothing: {fn}'

    # (5) penalty penetration: the platen sits P/(4 kn) below the column top
    # (the deck starts gap0 already closed, so the extra travel is that minus
    # gap0).  Bounded by P/(4 kn) AND equal to the closed form to 1e-6 rel.
    d = [ops.nodeDisp(t, 3) - ops.nodeDisp(p, 3) for t, p in zip(top, platen)]
    dmax = max(abs(x) for x in d)
    assert dmax <= G2_P / (4.0 * G2_KN) + 1.0e-12, (drain, d)
    dref = G2_P / (4.0 * G2_KN) - G2_GAP0
    assert max(abs(x - dref) for x in d) <= 1.0e-6 * dref, (drain, d, dref)

    # snapshot the contact model's state
    p_c = {n: ops.nodeDisp(n, 4) for n in nid.values()}
    uz_c = {n: ops.nodeDisp(n, 3) for n in nid.values()}

    # (4) the TWIN: equalDOF on DOFs 1-3, no contact
    rc2, nid2, top2, platen2 = _g2_run_bonded(drain)
    assert rc2 == 0, f'bonded twin did not converge (drain={drain})'
    p_b = {n: ops.nodeDisp(n, 4) for n in nid2.values()}
    uz_b = {n: ops.nodeDisp(n, 3) for n in nid2.values()}

    dp = max(abs(p_c[n] - p_b[n]) for n in p_c)
    assert dp <= 1.0e-8, f'pore pressure differs by {dp:.3e} (drain={drain})'
    # the drained static path: p == 0 everywhere in BOTH models
    assert max(abs(v) for v in p_c.values()) <= 1.0e-8, p_c
    assert max(abs(v) for v in p_b.values()) <= 1.0e-8, p_b

    # the column's vertical field is the same: the penalty acts only on the
    # platen's own extra travel, not on the column
    ref = max(abs(v) for v in uz_b.values())
    du = max(abs(uz_c[n] - uz_b[n]) for n in uz_c)
    assert du <= 1.0e-6 * ref, (
        f'column u_z differs by {du:.3e} (ref {ref:.3e}, drain={drain})')


def _g2_run_contact_reversed(drain, kn=G2_KN, P=G2_P, gap0=G2_GAP0):
    """The SAME deck with the contact roles swapped: the ndf-3 platen quad is
    the master facet and the four ndf-4 column top nodes are the SLAVES.

    This is the leg that exercises ADR-96 D3's headline site -- the NTS slave
    guard, which was a FATAL `getNumberDOF() != 3` (LadrunoContactHandler.cpp
    :1825) and is now `< 3`.  The `top` leg above only reaches the master-side
    skips.  Outward is -z: the platen's outward normal points down at the
    column, whose allowed half-space is below it."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 4)
    nid, top = _g2_column(drain)
    # start just PENETRATED, as the forward deck does: with outward -z the
    # slave's gap is (x_s - x_m).n = (L - z_platen)(-1), negative only for a
    # platen BELOW the column top.
    platen = _g2_platen(COL_L - gap0, fix_lateral=True)
    ops.contactSurface(10, '-master', 4, *platen)   # ndf-3 master facet
    ops.contactSurface(20, '-slave', *top)          # ndf-4 slaves
    ops.contact(1, 10, 20, kn, 0.0, 0.0, '-outward', 0.0, 0.0, -1.0)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for t in platen:
        ops.load(t, 0.0, 0.0, -P / 4.0)
    ops.constraints('LadrunoContact')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.integrator('LoadControl', 1.0)
    ops.test('NormDispIncr', 1.0e-12, 50, 0)
    ops.algorithm('Newton')
    ops.analysis('Static')
    rc = ops.analyze(1)
    fn = [ops.ladrunoContactForce(t) for t in top]
    return rc, nid, top, platen, fn


@pytest.mark.parametrize('drain', ['top', 'base'])
def test_g2_ndf4_nodes_on_the_nts_slave_side(drain):
    """G2, roles swapped: the ndf-4 column top nodes are the NTS SLAVES.
    Same three mechanical claims plus the twin's pore-pressure identity."""
    rc, nid, top, platen, fn = _g2_run_contact_reversed(drain)
    assert rc == 0, f'reversed-role contact model did not converge ({drain})'
    tot = sum(fn)
    assert abs(tot - G2_P) <= 1.0e-8 * G2_P, (drain, tot, fn)
    assert min(fn) > 0.0, (drain, fn)
    d = [ops.nodeDisp(p, 3) - ops.nodeDisp(t, 3) for t, p in zip(top, platen)]
    assert max(abs(x) for x in d) <= G2_P / (4.0 * G2_KN) + 1.0e-12, (drain, d)

    p_c = {n: ops.nodeDisp(n, 4) for n in nid.values()}
    uz_c = {n: ops.nodeDisp(n, 3) for n in nid.values()}
    rc2, nid2, top2, platen2 = _g2_run_bonded(drain)
    assert rc2 == 0
    p_b = {n: ops.nodeDisp(n, 4) for n in nid2.values()}
    uz_b = {n: ops.nodeDisp(n, 3) for n in nid2.values()}
    dp = max(abs(p_c[n] - p_b[n]) for n in p_c)
    assert dp <= 1.0e-8, f'pore pressure differs by {dp:.3e} (drain={drain})'
    ref = max(abs(v) for v in uz_b.values())
    du = max(abs(uz_c[n] - uz_b[n]) for n in uz_c)
    assert du <= 1.0e-6 * ref, (
        f'column u_z differs by {du:.3e} (ref {ref:.3e}, drain={drain})')
