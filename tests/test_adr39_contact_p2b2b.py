"""ADR-39 P2b-2b — `-kn auto`: automatic penalty stiffness from the master element.

`contact ... auto` sizes the penalty per (slave, master-segment) pair from the OWNING
solid element's initial stiffness, reduced through the segment normal:

    kn = f_si * mean_over_seg_nodes( nᵀ K_block_node n ),   f_si = 0.10

the element-stiffness form of LS-DYNA 26.14a (f·K·A²/V). Sourced generically through
base-Element virtuals (getExternalNodes / getInitialStiff) so any 3-DOF/node solid
works. Validated oracle-first in contact_prototypes/proto_p2b2b_autokn.py (6/6).

Tests (the slave-on-deformable-brick fixture from P2b-2a, now with `auto`):
  1. auto converges + lands a usable kn (penetration a small fraction of element size);
  2. kn ∝ E  (double-then-triple E → penetration scales 1/E; formulation-independent);
  3. ABSOLUTE match: P/pen == the oracle kn = f_si·mean(nᵀKn) recomputed from the SAME
     std trilinear-hex stiffness the C++ reads (imported from the validated prototype);
  4. no owning solid for the master face → pair skipped (warn), contact inert.

Static equilibrium of the loaded slave: contact force kn·|gap| balances P, so
|gap| = P/kn EXACTLY (independent of brick compression) ⇒ kn_auto = P/pen.
"""
import pathlib
import sys

import pytest

from _testbed import ops

# reuse the VALIDATED auto-kn oracle (single source of truth for the formula) so the
# absolute-match test ties the compiled reduction directly to proto_p2b2b_autokn.py.
_PROTO = pathlib.Path(__file__).resolve().parents[1] / "Ladruno_implementation" / "contact_prototypes"
sys.path.insert(0, str(_PROTO))
from proto_p2b2b_autokn import auto_kn, hex_stiffness, top_face  # noqa: E402

pytestmark = [pytest.mark.zone_a]


def _fixed_brick(L, E, nu, ele_tag=1, mat_tag=1):
    """8-node std LadrunoBrick on [0,L]^3, bottom face fixed. Returns top-4 node tags."""
    pts = [(0, 0, 0), (L, 0, 0), (L, L, 0), (0, L, 0),
           (0, 0, L), (L, 0, L), (L, L, L), (0, L, L)]
    tags = []
    for i, (x, y, z) in enumerate(pts):
        t = 1 + i
        ops.node(t, float(x), float(y), float(z))
        tags.append(t)
    for t in tags[:4]:
        ops.fix(t, 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", mat_tag, E, nu)
    ops.element("LadrunoBrick", ele_tag, *tags, mat_tag)   # default formulation = std
    return tags[4:]


def _run_autokn(L, E, nu, P):
    """slave node pressed onto a single fixed-bottom deformable brick with `-kn auto`.
    returns (penetration, top-corner node tag). Asserts convergence."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    top = _fixed_brick(L, E, nu)
    ops.node(99, L / 2, L / 2, L - 1.0e-8)     # just penetrated into the top face
    ops.fix(99, 1, 1, 0)                        # free z
    ops.contactSurface(10, "-master", 4, *top)
    ops.contactSurface(20, "-slave", 99)
    ops.contact(1, 10, 20, "auto", "-outward", 0.0, 0.0, 1.0)   # <-- auto penalty
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(99, 0.0, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-10, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "`-kn auto` static contact did not converge"
    z_slave = (L - 1.0e-8) + ops.nodeDisp(99, 3)
    z_top = L + ops.nodeDisp(top[0], 3)
    return z_top - z_slave, top


def test_autokn_converges_and_reasonable():
    """auto sizes a usable penalty: converges, and under a REPRESENTATIVE contact
    pressure the penetration is a small fraction of the element size (not ~O(L) =
    no-contact/blow-up, not locked). Penalty contact is inherently soft — the
    penetration scales with the applied pressure/modulus ratio; auto-kn (f_si=0.10,
    ≈ LS-DYNA's default) is deliberately a working default, not a hard constraint
    (stiffer contact = explicit `-kn $val`). So load with a working ~0.1%-modulus
    pressure P = p*E*L^2, where the penetration fraction is a few %."""
    L, E, nu = 1.0, 2.0e4, 0.0
    P = 1.0e-3 * E * L * L          # 0.1%-of-modulus contact pressure over the face
    pen, _ = _run_autokn(L, E, nu, P)
    assert pen > 0.0, f"no penetration (contact inert?) pen={pen}"
    kn_auto = P / pen
    assert kn_auto > 0.0
    assert pen / L < 0.10, f"penetration {pen/L:.3e} of element size too large (lock/blowup?)"


def test_autokn_scales_with_modulus():
    """kn ∝ E ⇒ penetration ∝ 1/E (formulation-independent, non-circular physics)."""
    L, nu, P = 1.0, 0.0, 1.0e2
    pen1, _ = _run_autokn(L, 1.0e4, nu, P)
    pen3, _ = _run_autokn(L, 3.0e4, nu, P)
    ratio = pen1 / pen3
    assert abs(ratio - 3.0) / 3.0 < 0.01, f"pen(E)/pen(3E)={ratio:.4f}, expected 3 (kn∝E)"


def test_autokn_matches_oracle_formula():
    """ABSOLUTE: P/pen == f_si·mean(nᵀKn) recomputed from the same std-hex stiffness
    (imported from the validated oracle). Pins the formula + the f_si=0.10 constant."""
    L, E, nu, P = 1.0, 2.0e4, 0.0, 1.0e2
    K, X = hex_stiffness(L=L, E=E, nu=nu)
    kn_oracle = auto_kn(K, top_face(X), n=[0.0, 0.0, 1.0])   # f_si=0.10 default
    pen, _ = _run_autokn(L, E, nu, P)
    kn_auto = P / pen
    rel = abs(kn_auto - kn_oracle) / kn_oracle
    assert rel < 0.02, f"kn_auto={kn_auto:.4f} vs oracle={kn_oracle:.4f} (rel {rel:.3e})"


def test_autokn_no_owning_solid_skips():
    """A master face whose nodes belong to NO 3-DOF/node solid → auto cannot size the
    penalty → the pair is skipped (warn) and contact is inert. The slave then settles
    on its own spring alone (disp = -P/k_spring), proving the contact contributed nothing."""
    L, P, KSPRING = 1.0, 1.0e2, 1.0e3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # 4 master face nodes at z=L, FIXED, NOT part of any solid element
    top = [1, 2, 3, 4]
    pts = [(0, 0, L), (L, 0, L), (L, L, L), (0, L, L)]
    for t, (x, y, z) in zip(top, pts):
        ops.node(t, float(x), float(y), float(z))
        ops.fix(t, 1, 1, 1)
    # slave just below the face (penetrated), anchored to ground by a z-spring so the
    # system is solvable; loaded DOWN (toward more penetration).
    ops.node(99, L / 2, L / 2, L - 1.0e-8)
    ops.fix(99, 1, 1, 0)
    ops.node(98, L / 2, L / 2, L - 1.0e-8)
    ops.fix(98, 1, 1, 1)                        # ground anchor
    ops.uniaxialMaterial("Elastic", 1, KSPRING)
    ops.element("zeroLength", 1, 98, 99, "-mat", 1, "-dir", 3)
    ops.contactSurface(10, "-master", 4, *top)
    ops.contactSurface(20, "-slave", 99)
    ops.contact(1, 10, 20, "auto", "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(99, 0.0, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-10, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "spring-only model did not converge"
    dz = ops.nodeDisp(99, 3)
    # contact inert (skipped) ⇒ pure spring response, no contact resistance
    assert abs(dz - (-P / KSPRING)) / (P / KSPRING) < 1e-6, \
        f"slave dz={dz} != -P/k={-P/KSPRING} (contact did NOT skip?)"
