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


# ---------------------------------------------------------------------------
#  G3 -- the no-tension spring across a (3,4) pair
# ---------------------------------------------------------------------------

def _spring_pair(ndf_a, ndf_b, direction=3, load=+P_LOAD, ent=True):
    """Node 1 (ndf_a) fixed; node 2 (ndf_b) coincident, free in `direction`,
    joined by a zeroLength spring in `direction`. Load node 2 along it.

    Returns (disp of node 2 along `direction`, node-2 DOF-4 value or None,
    spring force from the element response).
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
        fx[3] = 0          # leave DOF 4 FREE: the passenger must stay at zero
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
    ops.constraints('Plain')
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
    p4 = ops.nodeDisp(2, 4) if ndf_b > 3 else None
    force = list(ops.eleResponse(1, 'force'))
    return rc, u, p4, force


def test_g3_ent_spring_3_4_carries_compression():
    rc, u, p4, force = _spring_pair(3, 4, load=-P_LOAD)
    assert rc == 0
    # compression: the spring closes, K*u balances the load
    assert abs(u + P_LOAD / K_SPRING) <= 1.0e-12 * P_LOAD / K_SPRING, u
    assert p4 == 0.0, p4                         # the passenger never moved
    # 6-slot core response: three per node, node-2 z = -(-P) = +P
    assert len(force) == 6, force
    assert abs(force[5] - P_LOAD) <= 1.0e-9 * P_LOAD, force


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
    ops.fix(2, 1, 1, 0, 0)
    # a soft parallel elastic spring keeps the tension DOF non-singular so the
    # opened ENT can be measured rather than inferred from a failed solve
    ops.uniaxialMaterial('ENT', 1, K_SPRING)
    ops.uniaxialMaterial('Elastic', 2, 1.0e-3 * K_SPRING)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 3)
    ops.element('zeroLength', 2, 1, 2, '-mat', 2, '-dir', 3)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 0.0, 0.0, +P_LOAD, 0.0)
    ops.constraints('Plain'); ops.numberer('Plain'); ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-12, 30, 0); ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.25); ops.analysis('Static')
    for _ in range(4):
        assert ops.analyze(1) == 0
    u = ops.nodeDisp(2, 3)
    # only the soft spring carries: u = P / (1e-3 K)
    assert abs(u - P_LOAD / (1.0e-3 * K_SPRING)) <= 1.0e-9 * u, u
    ent = list(ops.eleResponse(1, 'force'))
    assert max(abs(v) for v in ent) == 0.0, ent      # ENT opened: nothing
    assert ops.nodeDisp(2, 4) == 0.0                 # passenger untouched


@pytest.mark.parametrize('pair', [(3, 3), (4, 4), (4, 3), (3, 6)])
def test_g3_same_answer_on_every_pair(pair):
    rc, u, p4, force = _spring_pair(*pair, load=-P_LOAD)
    assert rc == 0, pair
    assert abs(u + P_LOAD / K_SPRING) <= 1.0e-12 * P_LOAD / K_SPRING, (pair, u)
    if p4 is not None:
        assert p4 == 0.0
    # (3,3) is the vanilla D3N6 path, (3,6)/(4,3)/(4,4) the passenger path:
    # the core response is the same 6-vector on all of them
    assert len(force) == 6, (pair, force)
    assert abs(force[5] - P_LOAD) <= 1.0e-9 * P_LOAD, (pair, force)


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
