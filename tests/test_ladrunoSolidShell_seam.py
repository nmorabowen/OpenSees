"""LadrunoSolidShell (ADR 66 G10, ELE_TAG 33020) — the 6-DOF/3-DOF seam recipe.

An ASDShellQ4 panel butted to a LadrunoSolidShell patch: the mixed-dimension
modeling pattern the solid-shell exists for (shell everywhere, solid-shell in
the punching/bearing region). G10's charge: a documented, VALIDATED connection
recipe with the moment-continuity error quantified.

THE RECIPE (validated here, elastic exactness grade):

    ops.LadrunoTie('-shellSolid', '-slaveNodes', n, *solidFaceNodes,
                   '-masterEdge', nseg, a1, b1, ..., '-thickness', t)
    # static:   constraints('Lagrange')          (cross-ndf rows are
    #                                             Transformation-INCOMPATIBLE)
    # explicit: constraints('LadrunoProjection') + system('Diagonal') +
    #           solid-shell '-lumped' (G9) + small nodal rotary mass on the
    #           tied masters' theta DOFs (rigidLink -beam precedent, ADR-64)

MEASURED on the #504-era binary (matching-EI members, nu=0, imposed far-edge
rotation = pure constant-moment state):
  moment continuity   station M == theta*E*I/L on BOTH sides of the seam to
                      machine precision; interface rotation == theta/2 to 1e-9;
                      GP sigma_xx == M*z/I at EVERY solid GP (0.0% envelope
                      deviation) — the plane-section tie transmits the
                      constant-moment state EXACTLY between the element types
                      (no seam boundary layer at all in this state class)
  membrane            constant sigma_xx crosses at 1e-14 grade
  explicit            dt_cr-neutral to 1e-9 (no penalty stiffness), tie
                      violation 1.7e-17 over 300 CDL steps, transverse load
                      crosses into the -lumped solid-shell patch

REJECTED CANDIDATES (the reason each is not the recipe):
  equalDOF / translations-only rows   a moment HINGE: the shell rides freely
      (theta_edge == theta), the solid stays unbent, station M == 0. Fails
      loudly in the gates, silently in a real model.
  rigidLink 'beam'                    cross-ndf DEAD END, and a TRAP: vanilla
      RigidBeam refuses a 3-DOF constrained node vs 6-DOF retained node with
      an opserr WARNING ONLY — no exception, the constraint is simply not
      added, and the model solves DISCONNECTED (see LEDGER_quirks).
  LadrunoTie -hermite                 shell-edge-to-shell-edge (ndf-6 both
      sides) vocabulary; the -shellSolid/-hermite combo is refused outright
      (covered in tests/test_ladrunoTie_shellsolid.py REFUSE battery).

The generic tie machinery (nested/non-nested grids, |d|=0 identity rows, arm
guard, momentum conservation, BLOCKER preconditions) is certified in
tests/test_ladrunoTie_shellsolid.py against SSPbrick; this file certifies the
LadrunoSolidShell integration + the moment-continuity numbers.

Plan/ADR: Ladruno_implementation/66_ladruno_solidshell_adr.md (gate G10).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, NU, T, RHO = 1000.0, 0.0, 0.1, 1.0
C = T / 2.0
THETA = 0.01
I_W = T ** 3 / 12.0                      # bending inertia per unit width
M_REF = THETA * E * I_W / 2.0            # theta = M*L/(E*I), L = 2, width 1
XS = [0.0, 0.25, 0.5, 0.75, 1.0]         # solid-shell patch x-grid
ZS = [-C, C]                             # ONE element through thickness

_SHELL = {201: (1.0, 0.0), 202: (1.0, 0.5), 203: (1.0, 1.0),
          211: (2.0, 0.0), 212: (2.0, 0.5), 213: (2.0, 1.0)}
EDGE = [201, 202, 203]
FAR = [211, 212, 213]
QUADS = [[201, 211, 212, 202], [202, 212, 213, 203]]
MASTER_EDGE = [2, 201, 202, 202, 203]    # -masterEdge <nseg=2> a1 b1 a2 b2


def _st(i, j, k):
    return 1000 + i * 100 + j * 10 + k


def _mesh(ny=4, tie=True, lumped=False, rotary=False):
    """Solid-shell patch x in [0,1] (len(XS)-1 x ny x 1, t=0.1) + 2 ASDShellQ4
    x in [1,2]; the x=1 face nodes are the tie slaves. ny=4 nests 2:1 in the
    shell's y-split (consistent collocation); ny=2 matches it node-for-node
    (the rigidLink candidate needs coincident columns)."""
    ys = [j / ny for j in range(ny + 1)]
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('ElasticIsotropic', 1, E, NU, RHO)
    for i, x in enumerate(XS):
        for j, y in enumerate(ys):
            for k, z in enumerate(ZS):
                ops.node(_st(i, j, k), x, y, z)
    eid = 1
    solids = []
    for i in range(len(XS) - 1):
        for j in range(ny):
            conn = [_st(i, j, 0), _st(i + 1, j, 0),
                    _st(i + 1, j + 1, 0), _st(i, j + 1, 0),
                    _st(i, j, 1), _st(i + 1, j, 1),
                    _st(i + 1, j + 1, 1), _st(i, j + 1, 1)]
            opts = ('-lumped',) if lumped else ()
            ops.element('LadrunoSolidShell', eid, *conn, 1, '-nz', 2, *opts)
            solids.append(eid)
            eid += 1
    ops.model('basic', '-ndm', 3, '-ndf', 6)      # builder-ndf gotcha (ADR-64)
    for t, (x, y) in _SHELL.items():
        ops.node(t, x, y, 0.0)
    ops.section('ElasticMembranePlateSection', 1, E, NU, T, RHO)
    for e, conn in enumerate(QUADS, start=100):
        ops.element('ASDShellQ4', e, *conn, 1)
    if rotary:
        # explicit-projector rule (ADR-64 measured correction): tied masters'
        # theta DOFs stay in the explicit equation set -> small rotary mass
        for t in _SHELL:
            ops.mass(t, 0.0, 0.0, 0.0, 1e-4, 1e-4, 1e-4)
    if tie:
        face = _face(ny)
        ops.LadrunoTie('-shellSolid', '-slaveNodes', len(face), *face,
                       '-masterEdge', *MASTER_EDGE, '-thickness', T)
    return ys, solids


def _face(ny):
    return [_st(len(XS) - 1, j, k) for j in range(ny + 1) for k in range(2)]


def _clamp_root(ys):
    for j in range(len(ys)):
        for k in range(2):
            ops.fix(_st(0, j, k), 1, 1, 1)


def _static_solve():
    ops.constraints('Lagrange')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-12, 30)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    assert ops.analyze(1) == 0


def _station_NM(solids, ys, i):
    """Internal-force resultants (N, M about z=0) on the x-station i node set."""
    zmap = {_st(i, j, k): ZS[k] for j in range(len(ys)) for k in range(2)}
    N = M = 0.0
    for e in solids:
        nds = ops.eleNodes(e)
        f = list(ops.eleForce(e))
        for a, nd in enumerate(nds):
            if nd in zmap:
                N -= f[3 * a]
                M -= f[3 * a] * zmap[nd]
    return N, M


def _gp_sxx_envelope_dev(solids):
    """Max relative deviation of the per-element |sigma_xx| GP envelope from
    the beam form M*z/I at the Gauss fiber z = C/sqrt(3)."""
    exp = M_REF / I_W * C / 3.0 ** 0.5
    dev = 0.0
    for e in solids:
        s = list(ops.eleResponse(e, 'stress'))
        mx = max(abs(s[6 * g]) for g in range(len(s) // 6))
        dev = max(dev, abs(mx - exp) / exp)
    return dev


def _bend(ys):
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    _clamp_root(ys)
    for t in FAR:
        ops.sp(t, 5, THETA)              # imposed far-edge rotation: constant M
    _static_solve()
    return sum(ops.nodeDisp(m, 5) for m in EDGE) / 3.0


def _face_grads(ys):
    out = []
    for j in range(len(ys)):
        top = ops.nodeDisp(_st(len(XS) - 1, j, 1), 1)
        bot = ops.nodeDisp(_st(len(XS) - 1, j, 0), 1)
        out.append((top - bot) / (2.0 * C))
    return out


# --------------------------------------------------------------------------
# MEMBRANE — constant stress crosses the seam (nu=0, matching EA)
# --------------------------------------------------------------------------
def test_membrane_patch_crosses_seam():
    DELTA = 1e-3
    ys, solids = _mesh()
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    _clamp_root(ys)
    for t in FAR:
        ops.sp(t, 1, DELTA)
    _static_solve()
    sig = E * 0.5 * DELTA
    for s in _face(len(ys) - 1):
        assert abs(ops.nodeDisp(s, 1) - 0.5 * DELTA) <= 1e-10 * DELTA
    for e in solids:
        s = list(ops.eleResponse(e, 'stress'))
        for g in range(len(s) // 6):
            assert abs(s[6 * g] - sig) <= 1e-10 * sig, \
                f"element {e} lost the constant-stress state across the seam"


# --------------------------------------------------------------------------
# HEADLINE — moment continuity, quantified: the constant-moment state crosses
# the 6/3-DOF seam EXACTLY (elastic exactness, not a percent band)
# --------------------------------------------------------------------------
def test_moment_crosses_seam_exactly():
    ys, solids = _mesh()
    th_edge = _bend(ys)
    # kinematic continuity: uniform-EI member -> theta(x) = THETA*x/L linear;
    # a compliant or misaligned seam shows up here first
    assert th_edge == pytest.approx(0.5 * THETA, rel=1e-8), \
        "interface rotation != theta/2: the seam added compliance"
    # the face carries the full plane-section couple, uniformly across y
    for g in _face_grads(ys):
        assert g == pytest.approx(th_edge, rel=1e-8), \
            "face through-thickness gradient != interface rotation"
    # moment continuity: the station resultant equals theta*E*I/L on BOTH
    # sides of the seam (root column and seam-adjacent column)
    for i, sgn in ((0, +1.0), (len(XS) - 1, -1.0)):
        N, M = _station_NM(solids, ys, i)
        assert abs(N) <= 1e-10 * M_REF / T
        assert sgn * M == pytest.approx(M_REF, rel=1e-8), \
            f"station {i}: moment did not cross the seam intact"
    # stress-level continuity: EVERY solid GP sits on the beam form M*z/I —
    # no seam boundary layer in the constant-moment state class
    assert _gp_sxx_envelope_dev(solids) < 1e-8, \
        "GP sigma_xx deviated from M*z/I near the seam"


# --------------------------------------------------------------------------
# REJECTED CANDIDATE 1 — translations-only rows are a moment hinge
# --------------------------------------------------------------------------
def test_translations_only_is_hinge():
    ys, solids = _mesh(tie=False)
    for j, y in enumerate(ys):
        if y <= 0.5:
            (mA, mB), xi = (201, 202), y / 0.5
        else:
            (mA, mB), xi = (202, 203), (y - 0.5) / 0.5
        for k in range(2):
            for d in (1, 2, 3):
                args = [_st(len(XS) - 1, j, k), d, 1.0]
                for m, N in ((mA, 1.0 - xi), (mB, xi)):
                    if abs(N) > 1e-12:
                        args += [m, d, -N]
                ops.equationConstraint(*args)
    th_edge = _bend(ys)
    assert abs(th_edge - THETA) < 0.05 * THETA, \
        "hinge must ride freely (theta_edge == THETA)"
    _, M = _station_NM(solids, ys, len(XS) - 1)
    assert abs(M) < 1e-9 * M_REF, "a translations-only seam must carry no moment"
    assert max(abs(g) for g in _face_grads(ys)) < 0.05 * THETA, \
        "hinge leaked a through-thickness couple"


# --------------------------------------------------------------------------
# REJECTED CANDIDATE 2 — rigidLink 'beam' is a cross-ndf dead end AND a trap:
# RigidBeam warns (opserr) on the 3-vs-6 DOF mismatch, adds NOTHING, and the
# model solves DISCONNECTED. Pin the trap so the failure mode stays loud here.
# --------------------------------------------------------------------------
def test_rigidlink_cross_ndf_silently_disconnects():
    ys, solids = _mesh(ny=2, tie=False)          # node-matched seam columns
    for m, y in ((201, 0.0), (202, 0.5), (203, 1.0)):
        j = ys.index(y)
        for k in range(2):
            ops.rigidLink('beam', m, _st(len(XS) - 1, j, k))
    th_edge = _bend(ys)
    # the "connection" did nothing: free hinge + zero moment in the solid
    assert abs(th_edge - THETA) < 0.05 * THETA
    _, M = _station_NM(solids, ys, len(XS) - 1)
    assert abs(M) < 1e-9 * M_REF, \
        "rigidLink unexpectedly carried moment — the cross-ndf refusal changed"
    # the RECIPE on the same node-matched mesh is still exact
    ys, solids = _mesh(ny=2, tie=True)
    th_edge = _bend(ys)
    assert th_edge == pytest.approx(0.5 * THETA, rel=1e-8)


# --------------------------------------------------------------------------
# EXPLICIT — the recipe under CDL + LadrunoProjection + '-lumped' (G9) mass:
# dt_cr-neutral, machine-zero tie violation, the moment-carrying load crosses
# --------------------------------------------------------------------------
def _plane_section_pred(y, z):
    if y <= 0.5:
        (mA, mB), xi = (201, 202), y / 0.5
    else:
        (mA, mB), xi = (202, 203), (y - 0.5) / 0.5
    pred = [0.0, 0.0, 0.0]
    for m, N in ((mA, 1.0 - xi), (mB, xi)):
        u = [ops.nodeDisp(m, d) for d in (1, 2, 3)]
        th = [ops.nodeDisp(m, d) for d in (4, 5, 6)]
        pred[0] += N * (u[0] + th[1] * z)
        pred[1] += N * (u[1] - th[0] * z)
        pred[2] += N * u[2]
    return pred


def _explicit_analysis(handler):
    ops.constraints(handler)
    ops.numberer('Plain')
    ops.system('Diagonal')
    ops.test('NormDispIncr', 1e-12, 10)
    ops.algorithm('Linear')
    ops.integrator('CentralDifferenceLadruno', '-cfl')
    ops.analysis('Transient')


def test_explicit_seam_recipe():
    ys, solids = _mesh(tie=True, lumped=True, rotary=True)
    _clamp_root(ys)
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for t in FAR:
        ops.load(t, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0)
    _explicit_analysis('LadrunoProjection')
    assert ops.analyze(1, 1e-5) == 0
    dt_tied = ops.criticalTimeStep()
    dt = 0.4 * dt_tied
    viol = 0.0
    for _ in range(300):
        assert ops.analyze(1, dt) == 0
        for j, y in enumerate(ys):
            for k, z in enumerate(ZS):
                p = _plane_section_pred(y, z)
                s = _st(len(XS) - 1, j, k)
                for d in (1, 2, 3):
                    viol = max(viol, abs(ops.nodeDisp(s, d) - p[d - 1]))
    assert viol < 1e-10, f"kinematic tie must hold to machine zero (got {viol})"
    assert abs(ops.nodeDisp(_st(2, 0, 1), 3)) > 1e-4, \
        "the transverse load never crossed into the solid-shell patch"

    ys, solids = _mesh(tie=False, lumped=True, rotary=True)
    _explicit_analysis('Plain')
    assert ops.analyze(1, 1e-5) == 0
    dt_free = ops.criticalTimeStep()
    assert dt_tied == pytest.approx(dt_free, rel=1e-9), \
        "the projection tie must not touch dt_cr (no penalty stiffness)"
