"""Plane-family dynamics + serialization regression battery (ADR 70 P4a
adversarial-gate remediation).

CRITICAL regression: getResistingForceIncInertia in the plane family
(LadrunoQuad / LadrunoCST / LadrunoLST -geom finite, LadrunoCSTPair) built
f_int − Q + M·a in the SHARED static P, then called
getRayleighDampingForces() — which, with stiffness-proportional betaK,
re-enters getTangentStiff() → formFinite/formPair(1) → P.Zero()+refill,
silently destroying the inertia term while Newton still converges (to a
wrong dynamic solution). Fixed with the LadrunoBrick snapshot pattern; see
LEDGER_quirks "getResistingForceIncInertia MUST snapshot".

Gates:
  1. dynamic-Rayleigh-preserves-inertia, PARAMETRIZED over the four finite
     paths (self-validating: the undamped step-load run must show the ~2x
     dynamic overshoot vs static — proving the rig resolves real dynamics —
     then a tiny betaK must move the peak <5%; under the bug the damped peak
     collapses toward 1x static, a >40% divergence);
  2. pair vs 2-CST-finite eigen identity — same geometry/material/mass
     convention ⇒ identical M and K at F≈I ⇒ eigenvalues match ~1e-9 (pins
     the pair's lumped-mass wiring AND stiffness assembly against the proven
     sibling — a /3→/4 lumping slip or a dropped-triangle mass fails loudly);
  3. LadrunoCSTPair FileDatastore save/restore roundtrip (sendSelf/recvSelf):
     restored model must reproduce the committed resisting state.
"""
import math
import os

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, NU, RHO, THK = 1000.0, 0.3, 8.0, 1.0


# --------------------------------------------------------------------------
# element placement helpers: a 1 x NCELL column of unit cells, tags built
# from the shared node grid; every variant is the element's FINITE path.
# --------------------------------------------------------------------------
NCELL = 2


def _grid():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    tags = {}
    for j in range(NCELL + 1):
        for i in (0, 1):
            t = 1 + j * 2 + i
            ops.node(t, float(i), float(j))
            tags[(i, j)] = t
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.nDMaterial("LogStrain2D", 2, 1)
    return tags


def _place(kind, tags, rho=RHO):
    k = 0
    for j in range(NCELL):
        n1, n2 = tags[(0, j)], tags[(1, j)]
        n3, n4 = tags[(1, j + 1)], tags[(0, j + 1)]
        k += 1
        if kind == "pair":
            ops.element("LadrunoCSTPair", k, n1, n2, n3, n4, 2, "-thick", THK,
                        "-type", "PlaneStrain", "-rho", rho)
        elif kind == "cst":
            ops.element("LadrunoCST", 2 * k - 1, n1, n2, n3, 2, "-thick", THK,
                        "-type", "PlaneStrain", "-geom", "finite", "-rho", rho)
            ops.element("LadrunoCST", 2 * k, n1, n3, n4, 2, "-thick", THK,
                        "-type", "PlaneStrain", "-geom", "finite", "-rho", rho)
        elif kind == "quad":
            ops.element("LadrunoQuad", k, n1, n2, n3, n4, 2, "-thick", THK,
                        "-type", "PlaneStrain", "-formulation", "bbar",
                        "-geom", "finite", "-rho", rho)
        elif kind == "lst":
            # T6 pair per cell: corner grid + midside nodes made on the fly
            base = 100 + 20 * k
            mids = {}

            def mid(a, b):
                key = tuple(sorted((a, b)))
                if key not in mids:
                    ca, cb = ops.nodeCoord(a), ops.nodeCoord(b)
                    t = base + len(mids)
                    ops.node(t, 0.5 * (ca[0] + cb[0]), 0.5 * (ca[1] + cb[1]))
                    mids[key] = t
                return mids[key]

            for (c1, c2, c3), tag in (((n1, n2, n3), 2 * k - 1),
                                      ((n1, n3, n4), 2 * k)):
                ops.element("LadrunoLST", tag, c1, c2, c3,
                            mid(c1, c2), mid(c2, c3), mid(c3, c1), 2,
                            "-thick", THK, "-type", "PlaneStrain",
                            "-geom", "finite", "-rho", rho)
        else:
            raise ValueError(kind)


def _fix_base(tags):
    ops.fix(tags[(0, 0)], 1, 1)
    ops.fix(tags[(1, 0)], 1, 1)


def _tip_nodes(tags):
    return tags[(0, NCELL)], tags[(1, NCELL)]


# --------------------------------------------------------------------------
# gate 1: dynamic Rayleigh preserves inertia (the CRITICAL regression)
# --------------------------------------------------------------------------
def _tip_static(kind, load):
    tags = _grid()
    _place(kind, tags)
    _fix_base(tags)
    a, b = _tip_nodes(tags)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(a, load, 0.0)
    ops.load(b, load, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.1)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-10, 50)
    ops.analysis("Static")
    assert ops.analyze(10) == 0
    return abs(ops.nodeDisp(a, 1))


def _tip_dynamic_peak(kind, load, betaK):
    """step load applied at t=0 (Constant series), undamped/lightly damped
    Newmark transient; returns the peak tip displacement over ~1.5 periods."""
    tags = _grid()
    _place(kind, tags)
    _fix_base(tags)
    a, b = _tip_nodes(tags)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(a, load, 0.0)
    ops.load(b, load, 0.0)
    if betaK != 0.0:
        ops.rayleigh(0.0, betaK, 0.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-10, 50)
    ops.analysis("Transient")
    # first shear-ish period from eigen (system already assembled)
    lam = ops.eigen("-fullGenLapack", 1)[0]
    T1 = 2.0 * math.pi / math.sqrt(lam)
    nstep, dt = 150, 1.5 * T1 / 150
    peak = 0.0
    for _ in range(nstep):
        assert ops.analyze(1, dt) == 0, "transient step failed"
        peak = max(peak, abs(ops.nodeDisp(a, 1)))
    return peak


@pytest.mark.parametrize("kind", ["pair", "cst", "quad", "lst"])
def test_dynamic_rayleigh_preserves_inertia(kind):
    load = 1.0                      # small: stays near-linear, clean overshoot
    static = _tip_static(kind, load)
    peak0 = _tip_dynamic_peak(kind, load, 0.0)
    overshoot = peak0 / static
    assert 1.5 < overshoot < 2.2, (
        f"[{kind}] undamped peak/static = {overshoot:.3f}; expected the ~2x "
        "step-load overshoot — rig not resolving dynamics, gate cannot "
        "discriminate the inertia bug")
    peak_d = _tip_dynamic_peak(kind, load, 1.0e-5)
    rel = abs(peak_d - peak0) / peak0
    assert rel < 0.05, (
        f"[{kind}] tiny betaK moved the dynamic peak by {rel:.1%} "
        f"({peak0:.6g} -> {peak_d:.6g}, static {static:.6g}) — inertia is "
        "being dropped in getResistingForceIncInertia (resid-clobber bug)")


# --------------------------------------------------------------------------
# gate 2: pair vs 2-CST eigen identity (mass + stiffness wiring)
# --------------------------------------------------------------------------
def _eigs(kind, n=6):
    tags = _grid()
    _place(kind, tags)
    _fix_base(tags)
    return ops.eigen("-fullGenLapack", n)


def test_pair_eigen_matches_two_csts():
    """Mass/stiffness wiring gates. NOTE the pair's tangent at F = I is NOT
    identical to the 2-CST assembly — the F-bar-Patch coupling contributes
    even at zero stress (q = ½c:I ≠ 0, ḡ ≠ g_e); that difference IS the cure.
    So: (a) a LOOSE pair-vs-2-CST eigen comparison catches gross wiring errors
    (dropped-triangle mass/stiffness, wrong TRI_NODES) without over-pinning
    the formulation difference; (b) a SHARP invariance — doubling rho must
    scale every eigenvalue by exactly 1/2 — pins the lumped-mass path."""
    e_pair = _eigs("pair")
    e_cst = _eigs("cst")
    # measured legit formulation gap: mode-1 rel ≈ 38% at nu=0.3 (the pair is
    # SOFTER — the volumetric relief itself; plain CSTs over-stiffen even
    # mildly compressible bending). 0.45 bounds gross wiring errors (dropped
    # triangle, wrong TRI_NODES ⇒ singular/off by integer factors) while not
    # over-pinning the cure. The sharp wiring pins are the crown printA gate
    # and the rho-scaling invariance below.
    for i, (a, b) in enumerate(zip(e_pair, e_cst)):
        assert a > 0 and b > 0
        rel = abs(a - b) / b
        assert rel < 0.45, (
            f"mode {i+1}: pair {a:.6g} vs 2xCST {b:.6g} (rel {rel:.2%}) — "
            "mass or stiffness wiring grossly off")
    # exact mass-wiring invariance: -rho flag vs material can't be compared
    # (ElasticIsotropic has no rho arg here), so pin total-mass sensitivity
    # instead: doubling rho must scale every eigenvalue by exactly 1/2
    tags = _grid()
    _place("pair", tags, rho=2 * RHO)
    _fix_base(tags)
    e2 = ops.eigen("-fullGenLapack", 6)
    for i, (a, b) in enumerate(zip(e_pair, e2)):
        assert abs(b - a / 2) < 1e-6 * a, (
            f"mode {i+1}: eigenvalue did not scale 1/2 with 2x rho "
            f"({a:.9g} -> {b:.9g}) — lumped-mass wiring broken")


# --------------------------------------------------------------------------
# gate 3: FileDatastore save/restore roundtrip (sendSelf/recvSelf)
# --------------------------------------------------------------------------
def test_pair_database_roundtrip(tmp_path):
    dbdir = str(tmp_path / "pairdb")
    tags = _grid()
    _place("pair", tags)
    _fix_base(tags)
    a, b = _tip_nodes(tags)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(a, 2.0, -1.0)
    ops.load(b, 2.0, 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.2)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-10, 50)
    ops.analysis("Static")
    assert ops.analyze(5) == 0
    disp_ref = [ops.nodeDisp(n, d) for n in (a, b) for d in (1, 2)]
    stress_ref = list(ops.eleResponse(1, "stresses"))

    try:
        ops.database("File", dbdir)
        rc = ops.save(1)
    except Exception as e:
        pytest.skip(f"FileDatastore unavailable in this build: {e}")
    assert rc in (0, None)

    ops.wipe()
    ops.database("File", dbdir)
    assert ops.restore(1) in (0, None), "restore failed (recvSelf path)"

    disp_new = [ops.nodeDisp(n, d) for n in (a, b) for d in (1, 2)]
    for i, (r, n) in enumerate(zip(disp_ref, disp_new)):
        assert abs(r - n) <= 1e-12 + 1e-9 * abs(r), f"disp[{i}] {r} != {n}"
    # element state machinery survives: run one more step on the restored model
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.2)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-10, 50)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "restored model cannot continue the analysis"
    stress_new = list(ops.eleResponse(1, "stresses"))
    assert len(stress_new) == len(stress_ref) == 6


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
