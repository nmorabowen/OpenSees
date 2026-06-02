"""
Spring/ZeroLength behavioural-claim verification battery.

Empirically pins the *surprising, load-bearing* runtime claims from the
zerolength_and_link_springs_guide.md that source-reading alone cannot settle.
The pure code-level facts (the `-doRayleigh` default constants, the
`setTrialStrain` signatures) are confirmed by direct citation and are NOT
re-tested here -- they are transcription, not inference.

The headline ZeroLength `-doRayleigh 0` => zero stiffness-Rayleigh damping is
already pinned in test_damping_channels.py and is not duplicated.

What IS tested here (all via free-vibration log-decrement with implicit Newmark
average-acceleration, so there is no algorithmic damping to contaminate zeta,
except the two static / quasi-static probes):

  1. ZeroLengthSection `-doRayleigh` defaults ON (the OPPOSITE of plain
     ZeroLength): a zeroLengthSection built with NO flag damps under domain
     betaK Rayleigh, and built with `-doRayleigh 0` does NOT. (guide 2.2 / 3.2)
  2. ZeroLength + `Viscous` acts as a dashpot WITHOUT `-doRayleigh` and without
     any domain Rayleigh -- the absorbing-boundary recipe. (guide 4.3 / 5.2)
  3. A `Viscous` material ALONE gives a singular static tangent (its getTangent
     is 0); a static solve fails, and adding a parallel Elastic spring fixes it.
     (guide 5.2)
  4. THE DOUBLE-NEGATIVE: a Parallel(Elastic+Viscous) material damps inside a
     plain ZeroLength (rate is passed) but is INERT inside a
     SectionAggregator -> ZeroLengthSection (no rate channel on either side),
     even with Rayleigh disabled. (guide 2.7 / 4.2)
  5. TwoNodeLink applies its material damping with NO `-doRayleigh` and no domain
     Rayleigh (material getDampTangent is always added). (guide 2.5)
  6. CoupledZeroLength couples two directions through the deformation RESULTANT:
     an ElasticPP material pushed at 45 deg yields at resultant Fy, so each
     component force is Fy/sqrt(2) -- whereas two independent ZeroLength springs
     each carry the full Fy. (guide 2.4)

Run:  py -3.12 -m pytest tests/test_spring_damping_claims.py -v
(the dist/bin opensees.pyd is built for CPython 3.12)
"""
import os
import sys
import math

# --- bootstrap the locally-built engine -----------------------------------
DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
if os.path.isdir(DIST):
    try:
        os.add_dll_directory(DIST)
    except (FileNotFoundError, OSError):
        pass
    if DIST not in sys.path:
        sys.path.insert(0, DIST)

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

import numpy as np


# --- helpers ----------------------------------------------------------------
def _zeta_from_decay(t, u, n_skip=1, n_use=5):
    """Damping ratio of a free-decay record by log decrement (positive peaks)."""
    u = np.asarray(u)
    peaks = []
    for i in range(1, len(u) - 1):
        if u[i] > u[i - 1] and u[i] >= u[i + 1] and u[i] > 0.0:
            peaks.append((t[i], u[i]))
    assert len(peaks) >= n_skip + n_use + 1, (
        f"not enough peaks ({len(peaks)}) to measure damping"
    )
    p0 = peaks[n_skip][1]
    pn = peaks[n_skip + n_use][1]
    delta = math.log(p0 / pn) / n_use
    return delta / math.sqrt(4.0 * math.pi ** 2 + delta ** 2)


def _run_free(build, dt, nsteps, d0=1.0e-3, dof=1, rayleigh_domain=None, ndf=2):
    """Free-vibration transient. `build()` must create node 1 (fixed) and
    node 2 (mass, free in `dof`) plus the spring element. Returns (t[], u[]).

    ndf=2 for plain ZeroLength/TwoNodeLink; ndf=3 is REQUIRED for
    ZeroLengthSection (it works only on full 3-dof (2D) / 6-dof (3D) nodes)."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', ndf)
    build()
    ops.constraints('Transformation'); ops.numberer('Plain')
    ops.system('FullGeneral'); ops.test('NormDispIncr', 1e-12, 20)
    ops.algorithm('Linear'); ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')
    if rayleigh_domain is not None:
        ops.rayleigh(*rayleigh_domain)
    ops.setNodeDisp(2, dof, d0, '-commit')
    t, u = [], []
    for _ in range(nsteps):
        ops.analyze(1, dt)
        t.append(ops.getTime()); u.append(ops.nodeDisp(2, dof))
    return np.array(t), np.array(u)


# Common SDOF: f = 1 Hz, w = 2*pi, T = 1 s
M = 1.0
W = 2.0 * math.pi
K = W * W * M                      # ~39.478
DT = 0.005                        # 200 steps / period
NSTEPS = 1600                     # 8 periods
ZTGT = 0.05
C_VISC = 2.0 * ZTGT * math.sqrt(K * M)   # c for target zeta via dashpot
BETAK = 2.0 * ZTGT / W                    # betaK for target zeta via Rayleigh


# === 1. ZeroLengthSection -doRayleigh defaults ON ==========================
def _build_zls(do_rayleigh=None):
    def build():
        # ndf=3 REQUIRED for ZeroLengthSection; only dir 1 (axial/P) is free
        ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1, 1)
        ops.node(2, 0.0, 0.0); ops.fix(2, 0, 1, 1)
        ops.mass(2, M, 0.0, 0.0)
        ops.uniaxialMaterial('Elastic', 1, K)
        ops.section('Aggregator', 10, 1, 'P')      # axial response on local x
        if do_rayleigh is None:                    # no flag -> default (claim: ON)
            ops.element('zeroLengthSection', 1, 1, 2, 10)
        else:
            ops.element('zeroLengthSection', 1, 1, 2, 10,
                        '-doRayleigh', do_rayleigh)
    return build


def test_zeroLengthSection_doRayleigh_default_ON():
    # default (no flag) + domain betaK -> SHOULD damp (default is ON)
    t, u = _run_free(_build_zls(do_rayleigh=None), DT, NSTEPS,
                     rayleigh_domain=(0.0, 0.0, BETAK, 0.0), ndf=3)
    zeta = _zeta_from_decay(t, u)
    assert math.isclose(zeta, ZTGT, rel_tol=0.15), (
        f"zeroLengthSection default should damp (default -doRayleigh ON): "
        f"zeta={zeta:.4f}, target={ZTGT:.4f}")


def test_zeroLengthSection_doRayleigh_explicit_off_kills_damping():
    t, u = _run_free(_build_zls(do_rayleigh=0), DT, NSTEPS,
                     rayleigh_domain=(0.0, 0.0, BETAK, 0.0), ndf=3)
    zeta = _zeta_from_decay(t, u)
    assert abs(zeta) < 0.005, (
        f"zeroLengthSection -doRayleigh 0 should give ~zero damping, "
        f"got zeta={zeta:.4f}")


# === 2. ZeroLength + Viscous dashpot with NO -doRayleigh, NO Rayleigh ======
def test_zeroLength_viscous_damps_without_flag():
    def build():
        ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
        ops.node(2, 0.0, 0.0); ops.fix(2, 0, 1)
        ops.mass(2, M, 0.0)
        ops.uniaxialMaterial('Elastic', 1, K)
        ops.uniaxialMaterial('Viscous', 2, C_VISC, 1.0)
        # both materials on dir 1 (parallel); NO -doRayleigh (defaults 0)
        ops.element('zeroLength', 1, 1, 2, '-mat', 1, 2, '-dir', 1, 1)
    t, u = _run_free(build, DT, NSTEPS)        # no domain rayleigh at all
    zeta = _zeta_from_decay(t, u)
    assert math.isclose(zeta, ZTGT, rel_tol=0.15), (
        f"ZeroLength+Viscous should damp via stress path w/o any flag: "
        f"zeta={zeta:.4f}, target={ZTGT:.4f}")


# === 3. Viscous alone -> singular static tangent ===========================
def _static_unit_load(build):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    build()
    ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0, 0.0)
    ops.constraints('Transformation'); ops.numberer('Plain')
    ops.system('FullGeneral'); ops.test('NormDispIncr', 1e-8, 20)
    ops.algorithm('Newton'); ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    return ops.analyze(1)


def test_viscous_alone_singular_static_tangent():
    def build():
        ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
        ops.node(2, 0.0, 0.0); ops.fix(2, 0, 1)
        ops.uniaxialMaterial('Viscous', 2, C_VISC, 1.0)   # getTangent == 0
        ops.element('zeroLength', 1, 1, 2, '-mat', 2, '-dir', 1)
    rc = _static_unit_load(build)
    assert rc != 0, ("Viscous-only ZeroLength should give a SINGULAR static "
                     f"tangent (analyze should fail), got rc={rc}")


def test_viscous_parallel_elastic_static_ok():
    def build():
        ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
        ops.node(2, 0.0, 0.0); ops.fix(2, 0, 1)
        ops.uniaxialMaterial('Elastic', 1, K)
        ops.uniaxialMaterial('Viscous', 2, C_VISC, 1.0)
        ops.element('zeroLength', 1, 1, 2, '-mat', 1, 2, '-dir', 1, 1)
    rc = _static_unit_load(build)
    assert rc == 0, ("Viscous paralleled with Elastic should solve statically, "
                     f"got rc={rc}")


# === 4. Double-negative: Parallel(Elastic+Viscous) inert in ZeroLengthSection
def test_double_negative_viscous_in_zeroLengthSection_is_inert():
    # (a) plain ZeroLength: Parallel material's Viscous part gets a rate -> DAMPS
    def build_zl():
        ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
        ops.node(2, 0.0, 0.0); ops.fix(2, 0, 1)
        ops.mass(2, M, 0.0)
        ops.uniaxialMaterial('Elastic', 1, K)
        ops.uniaxialMaterial('Viscous', 2, C_VISC, 1.0)
        ops.uniaxialMaterial('Parallel', 3, 1, 2)
        ops.element('zeroLength', 1, 1, 2, '-mat', 3, '-dir', 1,
                    '-doRayleigh', 0)               # kill Rayleigh entirely
    t, u = _run_free(build_zl, DT, NSTEPS)
    zeta_zl = _zeta_from_decay(t, u)

    # (b) ZeroLengthSection consuming the SAME Parallel material via Aggregator:
    #     no strain rate reaches it -> Viscous part INERT -> ~undamped
    def build_zls():
        ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1, 1)        # ndf=3 for ZLS
        ops.node(2, 0.0, 0.0); ops.fix(2, 0, 1, 1)
        ops.mass(2, M, 0.0, 0.0)
        ops.uniaxialMaterial('Elastic', 1, K)
        ops.uniaxialMaterial('Viscous', 2, C_VISC, 1.0)
        ops.uniaxialMaterial('Parallel', 3, 1, 2)
        ops.section('Aggregator', 10, 3, 'P')
        ops.element('zeroLengthSection', 1, 1, 2, 10, '-doRayleigh', 0)
    t, u = _run_free(build_zls, DT, NSTEPS, ndf=3)
    zeta_zls = _zeta_from_decay(t, u)

    assert zeta_zl > 0.02, (
        f"sanity: Parallel(Elastic+Viscous) in a plain ZeroLength must damp, "
        f"got zeta_zl={zeta_zl:.4f}")
    assert abs(zeta_zls) < 0.005, (
        f"DOUBLE-NEGATIVE: the same Viscous in ZeroLengthSection must be INERT "
        f"(no rate channel), got zeta_zls={zeta_zls:.4f}")


# === 5. TwoNodeLink applies material damping with NO -doRayleigh ===========
def test_twoNodeLink_material_damping_without_flag():
    def build():
        # finite length 1.0 so local-x = global-x; node2 free in dir 1 only
        ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
        ops.node(2, 1.0, 0.0); ops.fix(2, 0, 1)
        ops.mass(2, M, 0.0)
        ops.uniaxialMaterial('Elastic', 1, K)
        ops.uniaxialMaterial('Viscous', 2, C_VISC, 1.0)
        ops.element('twoNodeLink', 1, 1, 2, '-mat', 1, 2, '-dir', 1, 1)
    t, u = _run_free(build, DT, NSTEPS)        # no -doRayleigh, no domain rayleigh
    zeta = _zeta_from_decay(t, u)
    assert math.isclose(zeta, ZTGT, rel_tol=0.15), (
        f"TwoNodeLink should apply material damping with no flag: "
        f"zeta={zeta:.4f}, target={ZTGT:.4f}")


# === 6. CoupledZeroLength couples through the deformation resultant =========
def _impose_45deg_disp(build, d):
    """Impose (d, d) at node 2 and return the element's node-2 force components.

    Uses the Penalty handler (NOT Transformation): prescribing both DOFs of
    node 2 via `sp` plus a fully-fixed node 1 leaves ZERO free equations, which
    the Transformation handler cannot solve and terminates the process. Penalty
    retains the DOFs (penalised), so the solve is well-posed."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    build()
    ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1)
    ops.sp(2, 1, d); ops.sp(2, 2, d)
    ops.constraints('Penalty', 1e14, 1e14); ops.numberer('Plain')
    ops.system('FullGeneral'); ops.test('NormDispIncr', 1e-8, 50)
    ops.algorithm('Newton'); ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    rc = ops.analyze(1)
    f = ops.eleForce(1)                 # [F1x, F1y, F2x, F2y]
    return rc, f[2], f[3]               # node-2 force components


def test_coupledZeroLength_yields_on_resultant():
    E = K
    epsy = 0.01
    Fy = E * epsy
    d = epsy                       # resultant sqrt(2)*epsy > epsy -> yielded

    # coupled: one ElasticPP on the resultant of dirs 1 & 2
    def build_coupled():
        ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
        ops.node(2, 0.0, 0.0)
        ops.uniaxialMaterial('ElasticPP', 1, E, epsy)
        ops.element('CoupledZeroLength', 1, 1, 2, 1, 2, 1)
    rc_c, c1, c2 = _impose_45deg_disp(build_coupled, d)

    # independent: two separate ElasticPP springs, one per direction
    def build_indep():
        ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
        ops.node(2, 0.0, 0.0)
        ops.uniaxialMaterial('ElasticPP', 1, E, epsy)
        ops.element('zeroLength', 1, 1, 2, '-mat', 1, 1, '-dir', 1, 2)
    rc_i, i1, i2 = _impose_45deg_disp(build_indep, d)

    assert rc_c == 0 and rc_i == 0, f"static solves failed: {rc_c}, {rc_i}"
    # coupled: resultant force capped at Fy -> each component Fy/sqrt(2)
    assert math.isclose(abs(c1), Fy / math.sqrt(2.0), rel_tol=0.05), \
        f"coupled comp |c1|={abs(c1):.4f}, expected Fy/sqrt2={Fy/math.sqrt(2):.4f}"
    assert math.isclose(abs(c1), abs(c2), rel_tol=0.02), "coupled comps unequal"
    res_c = math.hypot(c1, c2)
    assert math.isclose(res_c, Fy, rel_tol=0.05), \
        f"coupled resultant={res_c:.4f}, expected Fy={Fy:.4f}"
    # independent: each spring fully yields -> each component Fy
    assert math.isclose(abs(i1), Fy, rel_tol=0.05), \
        f"independent comp |i1|={abs(i1):.4f}, expected Fy={Fy:.4f}"
    # and the contrast: independent component is ~sqrt(2)x the coupled one
    assert abs(i1) > abs(c1) * 1.3, \
        f"expected independent ({abs(i1):.4f}) >> coupled ({abs(c1):.4f})"
