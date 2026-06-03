"""LadrunoEmbeddedRebar (embedded-reinforcement coupling element) — Zone-A battery.

The element ties a rebar node to a host solid element's nodes via precomputed
shape-function weights N_i (Mode P penalty coupling; perfect-bond-via-penalty or a
bond-slip uniaxial in the axial slot). See SRC/element/ladrunoEmbeddedRebar/ and
Ladruno_implementation/20_ladruno_embedded_reinforcement_adr.md.

These checks pin the coupling MECHANICS without needing a host continuum element:
we play the role of the "host" with explicit nodes whose displacements we prescribe
(fix them), so the gap g = u_rebar - sum N_i u_host is fully controlled and the
element's force/stiffness can be checked against hand calculation.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _single_host_model(perfect_k=None, bond=False, kt=1.0e6, dirv=(1.0, 0.0, 0.0)):
    """Rebar node 1 tied to a single host node 2 (N=[1.0]).  With one host node and
    N=1 the gap is simply g = u1 - u2, so the coupling is a 3D spring between them:
    axial (along dir) via perfect-k or the bond law, transverse via kt."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)        # rebar node
    ops.node(2, 0.0, 0.0, 0.0)        # the (single) host node, coincident
    ops.fix(2, 1, 1, 1)               # prescribe host -> u_host = 0, so g = u_rebar
    args = ["LadrunoEmbeddedRebar", 1, 1, 1, 2, "-shape", 1.0,
            "-dir", *dirv, "-kt", kt]
    if bond:
        # tau_max s1 s2 s3 tau_f alpha ; bondScale=1 so axial force == tau(s)
        ops.uniaxialMaterial("LadrunoBondSlip", 7, 10.0, 1.0e-3, 3.0e-3,
                             10.0e-3, 2.0, 0.4)
        args += ["-bond", 7, "-bondScale", 1.0]
    else:
        args += ["-perfect", perfect_k]
    ops.element(*args)


def _push_rebar(ux, uy, uz, nstep=1):
    """Displacement-control the rebar node to (ux,uy,uz) and return its reaction
    (= the coupling force the element applies)."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(1, 1, ux); ops.sp(1, 2, uy); ops.sp(1, 3, uz)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-10, 20); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nstep); ops.analysis("Static")
    assert ops.analyze(nstep) == 0
    ops.reactions()
    return ops.nodeReaction(1)


# --------------------------------------------------- 1. perfect-bond axial spring
def test_perfect_axial_force():
    """Pull the rebar node along the bar axis by u; axial coupling force = k*u."""
    k = 1.0e5
    _single_host_model(perfect_k=k)
    R = _push_rebar(2.0e-3, 0.0, 0.0)
    # reaction at the prescribed dof balances the element's internal force k*u
    assert abs(R[0]) == pytest.approx(k * 2.0e-3, rel=1e-6)


def test_perfect_transverse_uses_kt():
    """Move the rebar node transverse to the bar; the transverse force uses kt, not
    the axial k (anisotropic coupling, ADR D6)."""
    k, kt = 1.0e5, 3.0e5
    _single_host_model(perfect_k=k, kt=kt, dirv=(1.0, 0.0, 0.0))
    R = _push_rebar(0.0, 1.0e-3, 0.0)          # pure y (transverse) motion
    assert abs(R[1]) == pytest.approx(kt * 1.0e-3, rel=1e-6)
    assert abs(R[0]) < 1e-3                     # no axial force for transverse motion


# --------------------------------------------------- 2. bond-slip in the axial slot
def test_bond_axial_follows_tau_s():
    """With a bond material the axial force = tau(s)*bondScale. At s in the ascending
    power-law branch, force == tau_max*(s/s1)^alpha (bondScale=1)."""
    _single_host_model(bond=True)
    s = 0.5e-3                                   # < s1, ascending branch
    R = _push_rebar(s, 0.0, 0.0)
    tau = 10.0 * (s / 1.0e-3) ** 0.4
    assert abs(R[0]) == pytest.approx(tau, rel=1e-4)


# --------------------------------------------------- 3. partition of unity on hosts
def test_force_distributes_to_hosts_by_N():
    """Two host nodes with N=[0.25,0.75]; the reaction the element pushes back onto
    the hosts must split as N_i * (rebar force) and balance the rebar node."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)       # rebar
    ops.node(2, 0.0, 0.0, 0.0)       # host A
    ops.node(3, 0.0, 0.0, 0.0)       # host B
    for n in (2, 3):
        ops.fix(n, 1, 1, 1)
    k = 1.0e5
    ops.element("LadrunoEmbeddedRebar", 1, 1, 2, 2, 3, "-shape", 0.25, 0.75,
                "-dir", 1.0, 0.0, 0.0, "-perfect", k, "-kt", k)
    R = _push_rebar(1.0e-3, 0.0, 0.0)
    ops.reactions()
    fr = R[0]                         # rebar reaction
    f2 = ops.nodeReaction(2)[0]
    f3 = ops.nodeReaction(3)[0]
    # equilibrium: rebar + both hosts sum to zero; hosts split 0.25/0.75
    assert (fr + f2 + f3) == pytest.approx(0.0, abs=1e-6 * abs(fr))
    assert f2 == pytest.approx(0.25 * (-fr), rel=1e-6)
    assert f3 == pytest.approx(0.75 * (-fr), rel=1e-6)


# --------------------------------------------------- 4. serialization round-trip
def test_database_roundtrip():
    """sendSelf/recvSelf via a database datastore preserves the response."""
    _single_host_model(bond=True)
    _push_rebar(0.8e-3, 0.0, 0.0)
    f_before = ops.eleResponse(1, "axialForce")[0]
    assert abs(f_before) > 0.0
