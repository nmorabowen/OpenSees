"""LadrunoSolidShell (ADR 66 G7, ELE_TAG 33020) — moment-curvature vs LayeredShellFiberSection.

THE gate that any flexural claim must pass (ADR 66 D1 / ADR 19 D1): the same RC
wall section built two ways — the director-shell reference stack (ASDShellQ4 +
LayeredShell of LadrunoRCConcrete PlateFiber layers + PlateRebar steel) vs a
LadrunoSolidShell through-thickness stack (3D LadrunoRCConcrete at the GPs,
-nz lobatto thickness resolution, Steel01 trusses on the z=+-80 interface node
lines) — driven through the same section-level histories, with the deviation
QUANTIFIED. Section: t=200, strip b=100, rebar at z=+-80 (a_s=0.6 mm2/mm/face,
rho_tot=0.6%), fc=30/ft=3 harness backbones, cylindrical bending (eps_yy=0).

Two single-element section rigs (all DOF-prescribed / dummy-node driven, so the
comparison is section-true and solver-trivial):
  G7a  prescribed (eps0=0, kappa) ramp -> compare resultants N(k), M(k)
       [restrained wall bending; N develops as cracking is restrained]
  G7b  (u0, theta) dummy node + equationConstraint plane-section rows,
       N force-controlled, theta ramped -> the CLASSIC M-kappa at N=0 and at
       a wall service axial N = -0.1 fc Ag  [eps0 floats; N held exactly]

MEASURED on the #504-era binary (implicit both arms, kappa -> 6e-5 = deep
cracking + steel yield, 120 steps; deviations of solid vs shell):
  elastic anchor      both arms == E t^3/12(1-nu^2) closed form (solid EXACT;
                      the layered section is EXACT on its own midpoint-rule
                      quadrature, deficit Sum h_i^3/12 ~ 1-2% vs continuum)
  G7a stack6(nz 3)    dM = 2.1% of Mmax   dN = 16.4% of ft*t*b
  G7a 3-slice(nz 2/5/2, single core element) dN ~ 62% — the through-thickness
                      sigma_zz relaxation of ONE core element cannot follow the
                      cracked state; STACKING converges it (ADR 66 O2 answered:
                      element stacking, not per-GP-layer assignment)
  G7b N=0             dM = 2.2% of Mmax   de0 = 1.9%   N held to machine zero
  G7b N=-0.1 fc Ag    dM = 1.6% of Mmax   de0 = 5.1%   N held to machine zero

IMPLEX reporting quirk (pinned by test_implex_reporting_paths): under a fully
prescribed rig an -implex LadrunoRCConcrete inside ASDShellQ4 reports BIT-
IDENTICAL response to the implicit run — ASDShellQ4 output reflects the
post-commit state, and IMPL-EX re-integrates implicitly at commit; hosts that
report the trial state (LadrunoSolidShell, ShellMITC4) show the extrapolated
stresses (~3.5% here). On free-DOF problems implex DOES alter the ASDShellQ4
equilibrium path (0.7% measured on the G7b rig) — the flag is not inert, it is
invisible to prescribed-rig probes. The benchmark runs implicit on BOTH arms.

Plan/ADR: Ladruno_implementation/66_ladruno_solidshell_adr.md (gate G7, O2).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# ---- section constants (N, mm, MPa) ---------------------------------------
E, NU = 30000.0, 0.2
T, B, DX = 200.0, 100.0, 100.0          # thickness, strip width, rig length
ZR = 80.0                               # rebar plane |z| (20 mm cover to bar center)
AS_FACE = 0.6                           # steel area per mm width, per face
FY, ES, BS = 420.0, 200000.0, 0.01
FC, FT = 30.0, 3.0

# the rc_wall_harness backbones (strain-based; no -autoRegularization so both
# arms see the identical constitutive law regardless of element lch)
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0, 30.0, 5.0]
CD = [0.0, 0.0, 0.25, 1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010]
TS = [0.0, 3.0, 0.5]
TD = [0.0, 0.0, 1.0 - 0.5 / 5.0]

D_CYL = E * T ** 3 / 12.0 / (1.0 - NU ** 2)       # cylindrical plate rigidity / width

# solid stacks: interface planes + per-slice -nz (>=3 -> lobatto face points)
STACK6 = ([-100.0, -80.0, -40.0, 0.0, 40.0, 80.0, 100.0], [3, 3, 3, 3, 3, 3])
SLICE3 = ([-100.0, -80.0, 80.0, 100.0], [2, 5, 2])

KT, NSTEP = 6.0e-5, 120                 # ramp target (deep cracking + yield) / steps
DUM = 999


def _rc_concrete(tag, implex=False):
    args = ['LadrunoRCConcrete', tag, E, NU,
            '-Ce', *CE, '-Cs', *CS, '-Cd', *CD,
            '-Te', *TE, '-Ts', *TS, '-Td', *TD]
    if implex:
        args += ['-implex']
    ops.nDMaterial(*args)


# ---------------------------------------------------------------------------
# solid arm: 1x1 in-plane x n-slice LadrunoSolidShell stack + interface trusses
# ---------------------------------------------------------------------------
def _ntag(i, j, k):
    return 100 + i * 40 + j * 20 + k    # i x-station 0/1, j y 0/1, k z-level


def build_solid(stack=SLICE3, elastic=False, rebar=True, implex=False):
    zl, nzs = stack
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    if elastic:
        ops.nDMaterial('ElasticIsotropic', 1, E, NU)
    else:
        _rc_concrete(1, implex=implex)
    for i in range(2):
        for j in range(2):
            for k in range(len(zl)):
                ops.node(_ntag(i, j, k), i * DX, j * B, zl[k])
    eles = []
    for k in range(len(zl) - 1):
        conn = [_ntag(0, 0, k), _ntag(1, 0, k), _ntag(1, 1, k), _ntag(0, 1, k),
                _ntag(0, 0, k + 1), _ntag(1, 0, k + 1),
                _ntag(1, 1, k + 1), _ntag(0, 1, k + 1)]
        opts = ('-nz', nzs[k], '-quadz', 'lobatto') if nzs[k] >= 3 else ('-nz', nzs[k])
        ops.element('LadrunoSolidShell', k + 1, *conn, 1, *opts)
        eles.append(k + 1)
    if rebar:
        ops.uniaxialMaterial('Steel01', 11, FY, ES, BS)
        etag = 50
        for j in range(2):
            for zr in (-ZR, ZR):
                k = zl.index(zr)        # rebar planes must be stack interfaces
                ops.element('Truss', etag, _ntag(0, j, k), _ntag(1, j, k),
                            AS_FACE * B / 2.0, 11)
                eles.append(etag)
                etag += 1
    # cylindrical bending: uy = 0 everywhere; uz datum on one node
    for i in range(2):
        for j in range(2):
            for k in range(len(zl)):
                ops.fix(_ntag(i, j, k), 0, 1, 1 if (i, j, k) == (0, 0, 0) else 0)
    return eles, zl


def solid_resultants(eles, zl):
    """Section (N, M) about z=0 from internal forces on the x=0 station set."""
    zmap = {_ntag(0, j, k): zl[k] for j in range(2) for k in range(len(zl))}
    N = M = 0.0
    for e in eles:
        nds = ops.eleNodes(e)
        f = list(ops.eleForce(e))
        for a, nd in enumerate(nds):
            if nd in zmap:
                N -= f[3 * a]
                M -= f[3 * a] * zmap[nd]
    return N, M


# ---------------------------------------------------------------------------
# shell arm: 1 ASDShellQ4 on a LayeredShell section (RC) / elastic plate
# ---------------------------------------------------------------------------
SHN = {(0, 0): 1, (1, 0): 2, (1, 1): 3, (0, 1): 4}


def _shell_layers(ncov, ncore, rebar):
    """Layer stack [(matTag, thk), ...] bottom->top; rebar centroids at z=+-ZR."""
    tcov = T / 2.0 - ZR - (AS_FACE / 2.0 if rebar else 0.0)
    tcore = 2.0 * ZR - (AS_FACE if rebar else 0.0)
    cov = [(1, tcov / ncov)] * ncov
    core = [(1, tcore / ncore)] * ncore
    if rebar:
        return cov + [(21, AS_FACE)] + core + [(21, AS_FACE)] + cov
    return cov + core + cov


def build_shell(ncov=3, ncore=10, elastic=False, rebar=True, implex=False,
                elastic_rebar=False):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    for (i, j), n in SHN.items():
        ops.node(n, i * DX, j * B, 0.0)
    if elastic and not elastic_rebar:
        ops.section('ElasticMembranePlateSection', 10, E, NU, T)
        layers = None
    else:
        if elastic_rebar:
            ops.nDMaterial('ElasticIsotropic', 1, E, NU)
        else:
            _rc_concrete(1, implex=implex)
        if rebar:
            ops.uniaxialMaterial('Steel01', 11, FY, ES, BS)
            ops.nDMaterial('PlateRebar', 21, 11, 0.0)     # along the strip axis x
        layers = _shell_layers(ncov, ncore, rebar)
        flat = [v for lay in layers for v in lay]
        ops.section('LayeredShell', 10, len(layers), *flat)
    ops.element('ASDShellQ4', 1, 1, 2, 3, 4, 10)
    return layers


def shell_resultants():
    """Section (N, M) about z=0 from internal forces on the x=0 edge."""
    nds = ops.eleNodes(1)
    pos = {n: ij for ij, n in SHN.items()}
    N = M = 0.0
    f = list(ops.eleForce(1))
    for a, nd in enumerate(nds):
        if pos[nd][0] == 0:
            N -= f[6 * a]
            M -= f[6 * a + 4]
    return N, M


# ---------------------------------------------------------------------------
# solvers: G7a fully-prescribed ramp / G7b dummy-node section driver
# ---------------------------------------------------------------------------
def _analysis(handler, nsteps):
    ops.wipeAnalysis()
    if handler == 'Lagrange':
        ops.constraints('Lagrange')
        ops.numberer('Plain')
    else:
        ops.constraints('Transformation')
        ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-8, 50, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / nsteps)
    ops.analysis('Static')


def _cut_step(dlam, max_cuts=5):
    """One dlam advance with algorithm fallbacks + recursive halving."""
    ops.integrator('LoadControl', dlam)
    if ops.analyze(1) == 0:
        return True
    for alg in ('KrylovNewton', 'NewtonLineSearch'):
        ops.algorithm(alg)
        ok = ops.analyze(1)
        ops.algorithm('Newton')
        if ok == 0:
            return True
    if max_cuts == 0:
        return False
    return all(_cut_step(dlam / 2.0, max_cuts - 1) for _ in range(2))


def run_g7a_solid(kt=KT, nsteps=NSTEP, **kw):
    eles, zl = build_solid(**kw)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for i in range(2):
        for j in range(2):
            for k in range(len(zl)):
                ops.sp(_ntag(i, j, k), 1, kt * zl[k] * (i * DX - DX / 2.0))
    _analysis('Transformation', nsteps)
    out = []
    for _ in range(nsteps):
        if not _cut_step(1.0 / nsteps):
            break
        N, M = solid_resultants(eles, zl)
        out.append((ops.getTime() * kt, N, M))
    return out


def run_g7a_shell(kt=KT, nsteps=NSTEP, **kw):
    build_shell(**kw)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for (i, j), n in SHN.items():
        ops.sp(n, 1, 0.0)
        ops.sp(n, 5, kt * (i * DX - DX / 2.0))
        ops.fix(n, 0, 1, 1 if n == 1 else 0, 1, 0, 1)   # w free (datum on node 1)
    _analysis('Transformation', nsteps)
    out = []
    for _ in range(nsteps):
        if not _cut_step(1.0 / nsteps):
            break
        N, M = shell_resultants()
        out.append((ops.getTime() * kt, N, M))
    return out


def _staged_theta_ramp(kt, nsteps, n0, theta_cb, resultant_cb):
    """Phase A: ramp the axial preload (10 steps, then loadConst). Phase B:
    ramp the end rotation theta = kt*DX. Returns per-step resultant tuples."""
    if n0:
        _analysis('Lagrange', 10)
        for _ in range(10):
            assert ops.analyze(1) == 0, "axial preload phase failed"
        ops.loadConst('-time', 0.0)
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    theta_cb()
    _analysis('Lagrange', nsteps)
    out = []
    for _ in range(nsteps):
        if not _cut_step(1.0 / nsteps):
            break
        out.append(resultant_cb())
    return out


def run_g7b_solid(kt=KT, nsteps=NSTEP, n0=0.0, **kw):
    """(u0, theta) dummy-node plane-section driver: left station ux=0, right
    station tied to ux = u0 + theta*z via equationConstraint rows; axial force
    applied on u0 (work-conjugate), theta ramped by sp. eps0 floats."""
    eles, zl = build_solid(**kw)
    ops.node(DUM, DX, 0.0, 0.0)
    ops.fix(DUM, 0, 0, 1)
    for j in range(2):
        for k in range(len(zl)):
            ops.fix(_ntag(0, j, k), 1, 0, 0)
            args = [_ntag(1, j, k), 1, 1.0, DUM, 1, -1.0]
            if abs(zl[k]) > 1e-12:      # a zero rcoef is refused by the parser
                args += [DUM, 2, -zl[k]]
            ops.equationConstraint(*args)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    if n0:
        ops.load(DUM, n0, 0.0, 0.0)

    def theta():
        ops.sp(DUM, 2, kt * DX)

    def res():
        N, M = solid_resultants(eles, zl)
        return (ops.nodeDisp(DUM, 2) / DX, ops.nodeDisp(DUM, 1) / DX, N, M)

    return _staged_theta_ramp(kt, nsteps, n0, theta, res)


def run_g7b_shell(kt=KT, nsteps=NSTEP, n0=0.0, **kw):
    build_shell(**kw)
    for (i, j), n in SHN.items():
        if i == 0:
            ops.fix(n, 1, 1, 1, 1, 1, 1)              # clamped edge (theta ref)
        else:
            ops.fix(n, 0, 1, 0, 1, 0, 1)              # ux, w, theta_y free
    ops.equalDOF(2, 3, 1, 3, 5)                       # free edge moves as one
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    if n0:
        ops.load(2, n0, 0.0, 0.0, 0.0, 0.0, 0.0)

    def theta():
        ops.sp(2, 5, kt * DX)

    def res():
        N, M = shell_resultants()
        return (ops.nodeDisp(2, 5) / DX, ops.nodeDisp(2, 1) / DX, N, M)

    return _staged_theta_ramp(kt, nsteps, n0, theta, res)


# ---------------------------------------------------------------------------
# elastic anchor — both arms on the cylindrical closed form; the layered
# section's midpoint-rule quadrature is PREDICTED, not fudged
# ---------------------------------------------------------------------------
def test_g7_elastic_anchor():
    kt = 1.0e-5
    M_exact = D_CYL * kt * B

    out = run_g7a_solid(kt, 1, elastic=True, rebar=False)
    k, N, M = out[-1]
    assert M == pytest.approx(M_exact, rel=1e-6), \
        "solid stack must reproduce E t^3/12(1-nu^2) exactly (G2b closed form)"
    assert abs(N) < 1e-6 * M_exact / T

    out = run_g7a_shell(kt, 1, elastic=True, rebar=False)
    k, N, M = out[-1]
    assert M == pytest.approx(M_exact, rel=1e-6), \
        "ElasticMembranePlateSection reference must be exact too"

    # rebar bookkeeping, solid arm: concrete continuum + added steel, exact
    M_steel = ES * (AS_FACE * B) * ZR ** 2 * 2 * kt
    out = run_g7a_solid(kt, 1, elastic=True, rebar=True)
    assert out[-1][2] == pytest.approx(M_exact + M_steel, rel=1e-6), \
        "solid arm: trusses at z=+-80 add exactly Es*As*z^2 per face"

    # rebar bookkeeping, shell arm: EXACT on the layered section's own
    # midpoint-rule quadrature (one point per layer at its centroid) — the
    # deficit vs the continuum D is Sum h_i^3/12 per layer, NOT an error in
    # either element. This also pins the PlateRebar angle-0 = strip-axis
    # convention (a wrong angle would drop the steel term entirely).
    layers = build_shell(elastic=True, elastic_rebar=True, rebar=True)
    M_pred = 0.0
    zbot = -T / 2.0
    for mat, h in layers:
        zc = zbot + h / 2.0
        zbot += h
        Emod = ES if mat == 21 else E / (1.0 - NU ** 2)
        M_pred += Emod * h * B * zc ** 2 * kt
    out = run_g7a_shell(kt, 1, elastic=True, elastic_rebar=True, rebar=True)
    assert out[-1][2] == pytest.approx(M_pred, rel=1e-4), \
        "shell arm must sit exactly on its midpoint-rule layered prediction"
    assert 0.005 < 1.0 - M_pred / (M_exact + M_steel) < 0.03, \
        "the documented layered-quadrature deficit is ~1-2% at this layering"


# ---------------------------------------------------------------------------
# G7a — restrained (eps0 = 0) curvature ramp: resultant parity + the O2 answer
# ---------------------------------------------------------------------------
def test_g7a_restrained_section_parity():
    sh = run_g7a_shell(rebar=True)
    so = run_g7a_solid(stack=STACK6, rebar=True)
    assert len(sh) == NSTEP and len(so) == NSTEP, "both arms must complete"
    Mref = max(abs(m) for _, _, m in so)
    Nref = FT * T * B
    dM = max(abs(a[2] - b[2]) for a, b in zip(so, sh)) / Mref
    dN6 = max(abs(a[1] - b[1]) for a, b in zip(so, sh)) / Nref
    # measured 2.1% / 16.4% on the 6-element stack (see module docstring)
    assert dM < 0.04, f"restrained M(k) parity broke: {dM:.1%} of Mmax"
    assert dN6 < 0.25, f"restrained N(k) parity broke: {dN6:.1%} of ft*t*b"

    # O2 discriminator: ONE core element through [-80,80] cannot relax the
    # cracked-state sigma_zz profile — N(k) drifts; stacking converges it.
    so3 = run_g7a_solid(stack=SLICE3, rebar=True)
    assert len(so3) == NSTEP
    dM3 = max(abs(a[2] - b[2]) for a, b in zip(so3, sh)) / Mref
    dN3 = max(abs(a[1] - b[1]) for a, b in zip(so3, sh)) / Nref
    assert dM3 < 0.06, f"even the coarse stack holds M(k): {dM3:.1%}"
    assert dN3 > 1.5 * dN6, \
        "stacking must visibly converge the through-thickness state (O2: " \
        f"coarse dN {dN3:.1%} vs stack dN {dN6:.1%})"


# ---------------------------------------------------------------------------
# G7b — the classic moment-curvature (N held, eps0 free)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("n0, dm_tol, de_tol", [
    (0.0, 0.04, 0.05),                 # measured 2.2% / 1.9%
    (-0.1 * FC * T * B, 0.04, 0.10),   # wall axial ratio 0.1: measured 1.6% / 5.1%
], ids=["N0", "N_wall"])
def test_g7b_moment_curvature(n0, dm_tol, de_tol):
    so = run_g7b_solid(n0=n0, stack=STACK6, rebar=True)
    sh = run_g7b_shell(n0=n0, rebar=True)
    assert len(so) == NSTEP and len(sh) == NSTEP, "both arms must complete"
    # the rigs hold the axial resultant exactly (force-controlled on the
    # work-conjugate dummy/edge DOF) — this is what makes it an M-kappa rig
    assert max(abs(N - n0) for _, _, N, _ in so) < 1.0
    assert max(abs(N - n0) for _, _, N, _ in sh) < 1.0
    Mref = max(abs(m) for _, _, _, m in so)
    e0ref = max(abs(e) for _, e, _, _ in so) or 1.0
    dM = max(abs(a[3] - b[3]) for a, b in zip(so, sh)) / Mref
    de0 = max(abs(a[1] - b[1]) for a, b in zip(so, sh)) / e0ref
    assert dM < dm_tol, f"M-kappa parity broke at N0={n0:g}: {dM:.1%} of Mmax"
    assert de0 < de_tol, f"neutral-axis migration diverged: {de0:.1%}"
    # sanity: the ramp genuinely crossed cracking AND steel yield
    kappas = [k for k, _, _, _ in so]
    Ms = [m for _, _, _, m in so]
    EI0 = (Ms[0] / kappas[0])
    assert Ms[-1] / kappas[-1] < 0.35 * EI0, "section must be deeply cracked"


# ---------------------------------------------------------------------------
# IMPLEX reporting-path pin (see module docstring) — implicit-vs-implex on the
# PRESCRIBED rig: bit-identical under ASDShellQ4 (post-commit reporting; the
# IMPL-EX commit re-integrates implicitly), visibly different under the
# trial-state-reporting solid-shell. Neither is a defect; it means prescribed
# probes cannot see implex under ASDShellQ4 and parity gates must run
# implicit-vs-implicit (as the tests above do).
# ---------------------------------------------------------------------------
def test_implex_reporting_paths():
    kt, ns = 4.0e-5, 60
    a = run_g7a_shell(kt, ns, rebar=False, implex=False)
    b = run_g7a_shell(kt, ns, rebar=False, implex=True)
    assert len(a) == ns and len(b) == ns
    Mref = max(abs(m) for _, _, m in a)
    dsh = max(abs(x[2] - y[2]) for x, y in zip(a, b)) / Mref
    assert dsh < 1e-9, \
        f"ASDShellQ4 reports the post-commit (implicit) state: got rel {dsh:.2e}"
    a = run_g7a_solid(kt, ns, stack=SLICE3, rebar=False, implex=False)
    b = run_g7a_solid(kt, ns, stack=SLICE3, rebar=False, implex=True)
    assert len(a) == ns and len(b) == ns
    Mref = max(abs(m) for _, _, m in a)
    dso = max(abs(x[2] - y[2]) for x, y in zip(a, b)) / Mref
    assert dso > 1e-3, \
        "the solid-shell reports the trial state — implex must be visible " \
        f"(got rel {dso:.2e}); if this fails the reporting paths changed"
