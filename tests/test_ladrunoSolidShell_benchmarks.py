"""LadrunoSolidShell (ADR 66 G3, ELE_TAG 33020) — the classic shell obstacle course.

The two curved-geometry benchmarks every shell element must face, run on the
solid-shell (1 element through thickness, 2x2 x nz=2 gauss):

  G3a  Scordelis-Lo roof     cylindrical panel R=25, L=50, 40 deg, t=0.25,
                             E=4.32e8, nu=0, self-weight 90/area, rigid end
                             diaphragms; free-edge midpoint vertical deflection
                             reference 0.3024 (the standard FE normalization).
                             Membrane-dominated: the forgiving one.
  G3b  Pinched cylinder      R=300, L=600, t=3, E=3e6, nu=0.3, opposed point
                             loads P=1 at midspan, rigid end diaphragms;
                             radial deflection under the load reference
                             1.8248e-5. Inextensional bending + membrane
                             locking: the torture test.

These are the element's FIRST curved-mesh gates (P5.1 covered flat and
trapezoidal geometry): the zeta-fibers are radial, the faces are curved, and
the ANS tying operates off its exact flat-face class. Gates are (a) accuracy
bands at the canonical meshes set from measured values with margin, (b) a
monotone convergence trend under refinement, and (c) discrimination against
the '-formulation std' negative control (which locks on curved thin shells).

Both models exploit symmetry (quarter roof / one-eighth cylinder), all BCs
homogeneous (Transformation handler), single elastic LoadControl step.

Plan/ADR: Ladruno_implementation/66_ladruno_solidshell_adr.md (gate G3).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _model():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)


def _solve():
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-10, 20, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    assert ops.analyze(1) == 0


# ---------------------------------------------------------------------------
# G3a — Scordelis-Lo roof (quarter model)
# ---------------------------------------------------------------------------
_SL_R, _SL_L, _SL_PHI, _SL_T = 25.0, 50.0, math.radians(40.0), 0.25
_SL_E, _SL_NU, _SL_Q = 4.32e8, 0.0, 90.0          # q = self-weight per midsurface area
_SL_REF = 0.3024                                  # free-edge midpoint, downward


def _scordelis_ratio(n, formulation='ans'):
    """n x n quarter roof; returns |u_z(free-edge midpoint)| / 0.3024.

    Coordinates: theta from the crown (x = R sin th, z = R cos th), y axial.
    Quarter = theta in [0, 40deg], y in [0, L/2]. Symmetry: u_x = 0 at the
    crown plane (x=0), u_y = 0 at midspan (y=0); rigid diaphragm at y = L/2
    (u_x = u_z = 0); the theta = 40deg edge is free. Self-weight lumped per
    element midsurface area, split half to each face node."""
    _model()
    ops.nDMaterial('ElasticIsotropic', 1, _SL_E, _SL_NU)
    nt = ny = n
    dth, dy = _SL_PHI / nt, (_SL_L / 2.0) / ny
    tags = {}
    tag = 0
    for i in range(nt + 1):
        th = i * dth
        for j in range(ny + 1):
            for k, r in enumerate((_SL_R - _SL_T / 2, _SL_R + _SL_T / 2)):
                tag += 1
                tags[(i, j, k)] = tag
                ops.node(tag, r * math.sin(th), j * dy, r * math.cos(th))
    e = 0
    for i in range(nt):
        for j in range(ny):
            e += 1
            ops.element('LadrunoSolidShell', e,
                        tags[(i, j, 0)], tags[(i + 1, j, 0)],
                        tags[(i + 1, j + 1, 0)], tags[(i, j + 1, 0)],
                        tags[(i, j, 1)], tags[(i + 1, j, 1)],
                        tags[(i + 1, j + 1, 1)], tags[(i, j + 1, 1)],
                        1, '-formulation', formulation)
    for (i, j, k), t in tags.items():
        fx = 1 if i == 0 else 0                    # crown symmetry plane x = 0
        fy = 1 if j == 0 else 0                    # midspan symmetry plane y = 0
        fz = 0
        if j == ny:                                # rigid diaphragm (its own plane)
            fx, fz = 1, 1
        if fx or fy or fz:
            ops.fix(t, fx, fy, fz)
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    w_ele = _SL_Q * (_SL_R * dth * dy)             # element midsurface weight
    load = {}
    for i in range(nt):
        for j in range(ny):
            for (ii, jj) in ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)):
                for k in range(2):
                    t = tags[(ii, jj, k)]
                    load[t] = load.get(t, 0.0) + w_ele / 8.0
    for t, w in load.items():
        ops.load(t, 0.0, 0.0, -w)
    _solve()
    uz = 0.5 * (ops.nodeDisp(tags[(nt, 0, 0)])[2] + ops.nodeDisp(tags[(nt, 0, 1)])[2])
    return abs(uz) / _SL_REF


# ---------------------------------------------------------------------------
# G3b — pinched cylinder with end diaphragms (one-eighth model)
# ---------------------------------------------------------------------------
_PC_R, _PC_L, _PC_T = 300.0, 600.0, 3.0
_PC_E, _PC_NU, _PC_P = 3.0e6, 0.3, 1.0
_PC_REF = 1.8248e-5                                # radial, under the load


def _pinched_ratio(n, formulation='ans'):
    """n x n one-eighth cylinder; returns |u_y(under load)| / 1.8248e-5.

    Coordinates: x = R sin th, y = R cos th (theta from the loaded generator),
    z axial in [0, L/2]. Symmetry: u_x = 0 on the x = 0 plane (load plane),
    u_y = 0 on the y = 0 plane (theta = 90deg), u_z = 0 at midspan z = 0;
    rigid diaphragm at z = L/2 (u_x = u_y = 0). Load P/4 on the modeled
    eighth, split across the two face nodes, radially inward (-y)."""
    _model()
    ops.nDMaterial('ElasticIsotropic', 1, _PC_E, _PC_NU)
    nt = nz = n
    dth, dz = math.radians(90.0) / nt, (_PC_L / 2.0) / nz
    tags = {}
    tag = 0
    for i in range(nt + 1):
        th = i * dth
        for j in range(nz + 1):
            for k, r in enumerate((_PC_R - _PC_T / 2, _PC_R + _PC_T / 2)):
                tag += 1
                tags[(i, j, k)] = tag
                ops.node(tag, r * math.sin(th), r * math.cos(th), j * dz)
    e = 0
    for i in range(nt):
        for j in range(nz):
            e += 1
            # in-plane quad traversed (axial, then theta) so the bottom-face
            # normal points radially OUTWARD (detJ > 0 with k0 = inner face)
            ops.element('LadrunoSolidShell', e,
                        tags[(i, j, 0)], tags[(i, j + 1, 0)],
                        tags[(i + 1, j + 1, 0)], tags[(i + 1, j, 0)],
                        tags[(i, j, 1)], tags[(i, j + 1, 1)],
                        tags[(i + 1, j + 1, 1)], tags[(i + 1, j, 1)],
                        1, '-formulation', formulation)
    for (i, j, k), t in tags.items():
        fx = 1 if i == 0 else 0                    # load symmetry plane x = 0
        fy = 1 if i == nt else 0                   # theta = 90 symmetry plane y = 0
        fz = 1 if j == 0 else 0                    # midspan symmetry plane z = 0
        if j == nz:                                # rigid diaphragm (its own plane)
            fx, fy = 1, 1
        if fx or fy or fz:
            ops.fix(t, fx, fy, fz)
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for k in range(2):
        ops.load(tags[(0, 0, k)], 0.0, -_PC_P / 8.0, 0.0)   # P/4 on the eighth
    _solve()
    uy = 0.5 * (ops.nodeDisp(tags[(0, 0, 0)])[1] + ops.nodeDisp(tags[(0, 0, 1)])[1])
    return abs(uy) / _PC_REF


# ---------------------------------------------------------------------------
# gates — bands pinned to measured values (this binary, 2026-07-06) with
# margin:
#   Scordelis-Lo  ans: 0.9480 / 0.9739 / 0.9889 / 0.9923  (4/8/16/24)
#                 std 8x8: 0.1230 (locks by ~8x)
#   Pinched cyl   ans: 0.3781 / 0.7449 / 0.9271 / 0.9702  (4/8/16/24)
#                 std 8x8: 0.0693 (locks by ~11x)
# The pinched-cylinder 16x16 value 0.927 matches the canonical MITC4 result
# (~0.93 at 16x16) — the expected class behaviour for a transverse-shear/
# thickness-tied element WITHOUT in-plane membrane enhancement; the slow-from-
# below tail is the documented O3 (in-plane EAS) backlog signature, NOT a
# defect. Scordelis-Lo (membrane-dominated) sits ~1% off at 16x16 — inside
# the ADR's "~1-2% at the canonical meshes" claim.
# ---------------------------------------------------------------------------

def test_scordelis_lo_convergence():
    r = {n: _scordelis_ratio(n) for n in (4, 8, 16)}
    assert 0.96 < r[16] < 1.03, f"Scordelis-Lo 16x16 ratio {r[16]:.4f}"
    assert 0.94 < r[8] < 1.04, f"Scordelis-Lo 8x8 ratio {r[8]:.4f}"
    # refinement approaches the reference monotonically
    assert r[4] < r[8] < r[16], f"no monotone convergence trend: {r}"


def test_scordelis_lo_std_control():
    r_ans = _scordelis_ratio(8)
    r_std = _scordelis_ratio(8, formulation='std')
    assert r_std < 0.3 * r_ans, \
        f"discrimination lost: std {r_std:.4f} vs ans {r_ans:.4f}"


def test_pinched_cylinder_convergence():
    r = {n: _pinched_ratio(n) for n in (8, 16, 24)}
    assert 0.88 < r[16] < 1.02, f"pinched cylinder 16x16 ratio {r[16]:.4f}"
    assert 0.93 < r[24] < 1.02, f"pinched cylinder 24x24 ratio {r[24]:.4f}"
    assert r[8] < r[16] < r[24], \
        f"refinement must soften monotonically toward the reference: {r}"


def test_pinched_cylinder_std_control():
    r_ans = _pinched_ratio(8)
    r_std = _pinched_ratio(8, formulation='std')
    assert r_std < 0.2 * r_ans, \
        f"discrimination lost: std {r_std:.4f} vs ans {r_ans:.4f}"
