"""ADR-30 Phase 2 — general-C transport (rigidLink / rigidDiaphragm) Zone-A battery.

P1 proved the projector operator is general-C (L=[I;C]); the handler builds L from
the MP's full Ccr matrix. So rigidLink -beam and rigidDiaphragm transport constraints
flow through the SAME code with NO new projector logic — P2 VALIDATES that and pins the
boundary conditions (tied DOFs need lumped mass; tied DOFs must not be SP-fixed).

  T6fix  the P0/T6 transport model (offset slave via rigidLink -beam) under
         LadrunoProjection+Diagonal now MATCHES the dense-correct Transformation+
         FullGeneral reference — i.e. the projection handler FIXES the silent
         off-diagonal-mass drop P0 documented. Transport gap u_c-(C u_r) ~ machine eps.
  T7     rigidLink -beam free vibration: the projected (Diagonal) run matches the
         assembled (FullGeneral) reference frequency/response (the transport rotary
         inertia is carried correctly).
  T3     3D rigidDiaphragm slab + corner masses, lateral+torsional load: projection vs
         Transformation+FullGeneral; the diaphragm stays rigid to machine eps.
  BOUNDARY  a massless tied (rotational) slave DOF is refused (the §2.5 limit for
            rotational ties); SP-fixing a DOF the diaphragm controls is refused.

All run against the shipped P1 build — general C needs no rebuild.
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"
PROJ = "LadrunoProjection"


def _rms(a):
    return float(np.sqrt(np.mean(a * a)))


def _rel(a, b):
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    d = _rms(a) + _rms(b)
    return 0.0 if d == 0.0 else _rms(a - b) / (0.5 * d)


# --------------------------------------------------------------------------
# 2D offset-slave-via-rigidLink-beam model (the P0/T6 transport case)
# --------------------------------------------------------------------------
_M0, _IRZ, _MS, _MSRZ, _D = 1.0, 1.0, 1.0, 0.2, 2.0
_KX = _KY = _KR = 100.0


def _run_rigidlink_beam(handler, system, integrator, moment=5.0, nsteps=400, dt=0.01):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(10, 0.0, 0.0); ops.fix(10, 1, 1, 1)
    ops.node(1, 0.0, 0.0); ops.mass(1, _M0, _M0, _IRZ)
    ops.node(2, _D, 0.0);  ops.mass(2, _MS, _MS, _MSRZ)   # slave needs rotational mass
    for i, k in ((1, _KX), (2, _KY), (3, _KR)):
        ops.uniaxialMaterial("Elastic", i, k)
    ops.element("zeroLength", 1, 10, 1, "-mat", 1, 2, 3, "-dir", 1, 2, 3)
    ops.rigidLink("beam", 1, 2)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, moment)                          # moment on master rz, from rest
    ops.constraints(*handler); ops.numberer("Plain"); ops.system(*system)
    ops.test("NormDispIncr", 1e-10, 50)
    if integrator == "cdl":
        ops.algorithm("Linear"); ops.integrator(CDL)
    else:
        ops.algorithm("Newton"); ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    uy = [ops.nodeDisp(1, 2)]
    gap = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return None, None
        uy.append(ops.nodeDisp(1, 2))
        gap = max(gap, abs(ops.nodeDisp(2, 2) - (ops.nodeDisp(1, 2) + _D * ops.nodeDisp(1, 3))))
    return np.array(uy), gap


def test_T6fix_projection_fixes_the_transport_mass_drop():
    """The headline P2 result: LadrunoProjection+Diagonal reproduces the dense-correct
    Transformation+FullGeneral answer on the exact case where Transformation+Diagonal
    silently dropped the condensed off-diagonal mass (P0/T6)."""
    proj, gap = _run_rigidlink_beam((PROJ,), ("Diagonal",), "cdl")
    ref, _ = _run_rigidlink_beam(("Transformation",), ("FullGeneral",), "cdl")
    imp, _ = _run_rigidlink_beam(("Transformation",), ("FullGeneral",), "imp")
    assert proj is not None and ref is not None, "a run failed unexpectedly"
    assert _rms(proj) > 1e-3, "transport-induced uy too small to discriminate"
    assert _rel(proj, ref) < 1e-6, (
        f"projection does not match the dense reference (rel {_rel(proj, ref):.3e}) — "
        "general-C transport is wrong")
    assert _rel(proj, imp) < 0.1, f"projection vs implicit too far ({_rel(proj, imp):.3e})"
    assert gap < 1e-12, f"transport constraint u_c=(u_r+d*rz) not held ({gap:.3e})"


def _analytic_uy_p2_forced(t, moment=5.0):
    """OpenSees-INDEPENDENT golden truth (Gate-C): closed-form forced modal response of
    the condensed (uy,rz) master system for the rigidLink-beam model, hand-assembled from
    the physics (NOT via any OpenSees handler). The condensed mass includes the slave's
    offset translational mass (transport) AND the slave's own rz mass tied to the master:
        M22 = [[M0+MS,        MS*D            ],
               [MS*D,  IRZ + MS*D^2 + MSRZ    ]],   K22 = diag(KY, KR)
    Forced from rest by a constant moment on rz: u(t) = u_s + sum_i phi_i*eta_i(t),
    u_s = K22^-1 F, eta_i(t) = -(phi_i^T M u_s) cos(w_i t)."""
    t = np.asarray(t, dtype=float)
    M22 = np.array([[_M0 + _MS, _MS * _D],
                    [_MS * _D, _IRZ + _MS * _D * _D + _MSRZ]])
    K22 = np.array([[_KY, 0.0], [0.0, _KR]])
    F = np.array([0.0, moment])
    u_s = np.linalg.solve(K22, F)
    w2, phi = np.linalg.eig(np.linalg.solve(M22, K22))
    w = np.sqrt(w2.real); phi = phi.real
    for i in range(phi.shape[1]):
        phi[:, i] /= math.sqrt(phi[:, i] @ M22 @ phi[:, i])     # M-orthonormalize
    uy = np.full_like(t, u_s[0])
    for i in range(len(w)):
        eta0 = -(phi[:, i] @ M22 @ u_s)
        uy += phi[0, i] * eta0 * np.cos(w[i] * t)
    return uy


def test_T6b_projection_matches_analytic_independent_reference():
    """Gate-C rigor: close the validation triangle with an OpenSees-INDEPENDENT anchor.
    P0 showed Transformation+Full ~ analytic; T6fix shows proj ~ Transformation+Full;
    T6b shows proj ~ analytic directly, so the proj/Transformation agreement is not a
    shared-Ccr artifact. The analytic includes the slave rotational mass (transport)."""
    dt, nsteps = 0.01, 400
    proj, _ = _run_rigidlink_beam((PROJ,), ("Diagonal",), "cdl", moment=5.0, nsteps=nsteps, dt=dt)
    assert proj is not None
    anchor = _analytic_uy_p2_forced(np.arange(nsteps + 1) * dt, moment=5.0)
    assert _rel(proj, anchor) < 2e-2, (
        f"projection does not match the handler-independent analytic reference "
        f"(rel {_rel(proj, anchor):.3e}) — transport math is suspect")


def test_T7_rigidlink_beam_matches_assembled_reference():
    """rigidLink -beam free-vibration content: the Diagonal projected run matches the
    FullGeneral assembled run (the rotary transport inertia is carried right)."""
    # smaller moment + longer record to exercise the oscillation, not a ramp
    proj, _ = _run_rigidlink_beam((PROJ,), ("Diagonal",), "cdl", moment=2.0, nsteps=600)
    ref, _ = _run_rigidlink_beam(("Transformation",), ("FullGeneral",), "cdl",
                                 moment=2.0, nsteps=600)
    assert proj is not None and ref is not None
    assert _rel(proj, ref) < 1e-6, f"rel {_rel(proj, ref):.3e}"


# --------------------------------------------------------------------------
# T3 — 3D rigid diaphragm
# --------------------------------------------------------------------------
def _run_diaphragm(handler, system, integrator, nsteps=300, dt=0.005):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    # master at the centroid: in-plane translation + perpendicular (rz) rotation free
    ops.node(100, 0.0, 0.0, 0.0)
    ops.mass(100, 4.0, 4.0, 1e-9, 1e-9, 1e-9, 8.0)
    ops.fix(100, 0, 0, 1, 1, 1, 0)
    corners = [(1, 1.0, 1.0), (2, -1.0, 1.0), (3, -1.0, -1.0), (4, 1.0, -1.0)]
    gnd = 10
    for tag, x, y in corners:
        ops.node(tag, x, y, 0.0)
        # corner: ux,uy translational mass + a SMALL rz rotational mass (tied DOF needs
        # nonzero lumped mass under the projection handler — the §2.5 boundary)
        ops.mass(tag, 1.0, 1.0, 1e-9, 1e-9, 1e-9, 0.01)
        ops.fix(tag, 0, 0, 1, 1, 1, 0)                # free ux,uy,rz (rz is diaphragm-tied)
        ops.node(gnd, x, y, -1.0); ops.fix(gnd, 1, 1, 1, 1, 1, 1)
        ops.uniaxialMaterial("Elastic", tag, 50.0)
        ops.element("zeroLength", tag, gnd, tag, "-mat", tag, tag, "-dir", 1, 2)
        gnd += 1
    ops.rigidDiaphragm(3, 100, 1, 2, 3, 4)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
    ops.load(100, 3.0, 0.0, 0.0, 0.0, 0.0, 2.0)       # lateral + torque, from rest
    ops.constraints(*handler); ops.numberer("Plain"); ops.system(*system)
    ops.test("NormDispIncr", 1e-10, 50)
    if integrator == "cdl":
        ops.algorithm("Linear"); ops.integrator(CDL)
    else:
        ops.algorithm("Newton"); ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    ux = [ops.nodeDisp(100, 1)]
    gap = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return None, None
        ux.append(ops.nodeDisp(100, 1))
        # corner1 (x=1,y=1): rigid-diaphragm ux = master.ux - master.rz*(y1-ym)
        gap = max(gap, abs(ops.nodeDisp(1, 1) - (ops.nodeDisp(100, 1) - ops.nodeDisp(100, 6) * 1.0)))
    return np.array(ux), gap


def test_T3_rigid_diaphragm_matches_reference_and_stays_rigid():
    proj, gap = _run_diaphragm((PROJ,), ("Diagonal",), "cdl")
    ref, _ = _run_diaphragm(("Transformation",), ("FullGeneral",), "cdl")
    assert proj is not None, "projection diaphragm run failed"
    assert ref is not None, "reference diaphragm run failed"
    assert _rms(proj) > 1e-3
    assert _rel(proj, ref) < 1e-6, (
        f"diaphragm projection vs dense reference rel {_rel(proj, ref):.3e}")
    assert gap < 1e-10, f"diaphragm not rigid (gap {gap:.3e})"


# --------------------------------------------------------------------------
# BOUNDARY — the documented limits for transport ties
# --------------------------------------------------------------------------
def _runs_ok(build):
    try:
        build()
        return ops.analyze(1, 1e-3) == 0
    except Exception:
        return False


def test_massless_rotational_tie_refused():
    """A rigidLink -beam tying a slave rz that has ZERO rotational mass = a massless
    tied DOF -> refused (the §2.5 boundary; with Transformation the slave is eliminated
    so it is fine, but the projection keeps the slave's own equation)."""
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(10, 0.0, 0.0); ops.fix(10, 1, 1, 1)
        ops.node(1, 0.0, 0.0); ops.mass(1, 1.0, 1.0, 1.0)
        ops.node(2, 2.0, 0.0); ops.mass(2, 1.0, 1.0, 0.0)    # slave rz massless
        for i, k in ((1, 100.0), (2, 100.0), (3, 100.0)):
            ops.uniaxialMaterial("Elastic", i, k)
        ops.element("zeroLength", 1, 10, 1, "-mat", 1, 2, 3, "-dir", 1, 2, 3)
        ops.rigidLink("beam", 1, 2)
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(1, 0., 0., 5.)
        ops.constraints(PROJ); ops.numberer("Plain"); ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-10, 50); ops.algorithm("Linear")
        ops.integrator(CDL); ops.analysis("Transient")
    assert not _runs_ok(build), "massless rotational tie should be refused"


def test_sp_fixed_on_diaphragm_tied_dof_refused():
    """Fixing a DOF the diaphragm controls (slave rz) = SP-on-slave -> refused."""
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(100, 0., 0., 0.); ops.mass(100, 4., 4., 1e-9, 1e-9, 1e-9, 8.)
        ops.fix(100, 0, 0, 1, 1, 1, 0)
        ops.node(1, 1., 1., 0.); ops.mass(1, 1., 1., 1e-9, 1e-9, 1e-9, 0.01)
        ops.fix(1, 0, 0, 1, 1, 1, 1)                # rz FIXED but diaphragm ties it
        ops.node(10, 1., 1., -1.); ops.fix(10, 1, 1, 1, 1, 1, 1)
        ops.uniaxialMaterial("Elastic", 1, 50.)
        ops.element("zeroLength", 1, 10, 1, "-mat", 1, 1, "-dir", 1, 2)
        ops.rigidDiaphragm(3, 100, 1)
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(100, 3., 0., 0., 0., 0., 2.)
        ops.constraints(PROJ); ops.numberer("Plain"); ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-10, 50); ops.algorithm("Linear")
        ops.integrator(CDL); ops.analysis("Transient")
    assert not _runs_ok(build), "SP-fixed diaphragm-tied DOF should be refused"
