"""WP-160 -- ``ManzariDafalias`` ``MaxStrainInc`` / ``MaxEnergyInc`` sub-stepped
with UNINITIALISED moduli and never advanced the elastic strain.

THE DEFECTS (vanilla). Both sub-steppers declared ``double nDGamma,
nVoidRatio, nG, nK;`` and handed ``nG, nK`` by reference to every sub-step as
its elastic moduli; ``ForwardEuler`` and ``ModifiedEuler`` build ``aC =
GetStiffness(K, G)`` from them before writing them. And both loops advanced
``cStress, cStrain, cAlpha, cFabric`` but never ``cEStrain``, so each
sub-step restarted from the INITIAL elastic strain. Reach:

  * ``MaxStrainInc`` -- IntScheme 7, 8, 9 (its switch sends EVERY case to
    ``ForwardEuler``; the names MFE / RK are not honoured -- owner decision,
    WP-160: unchanged), whenever the largest strain component of the
    increment exceeds ``maxStrainInc = 1e-5``;
  * ``MaxEnergyInc`` -- IntScheme 4 (FE halves), 0 (ModifiedEuler halves),
    6 (RungeKutta4 halves), whenever ``de : dsigma > 1e-4``. RK4 WRITES the
    G, K it is handed, so the fix also keeps the entry moduli for the halves.

WP-160 gives every sub-step the moduli the un-sub-stepped call gets (the
committed G, K) and advances ``cEStrain``. Sub-stepping is then what it
claims to be: a refinement of the single-step call with the same integrator.

GATES (one increment from a committed plastic state ON the yield surface --
history under IntScheme 1 -- then the scheme under test for the last step
only, switched by ``setParameter IntegrationScheme``; built as IntScheme 1
throughout, so the constructor's scheme-3/5 warning latch is never touched):

  1. ``test_forward_euler_substeps_keep_the_elastic_bookkeeping`` [4, 7, 8,
     9] -- with the sub-step moduli frozen at the committed Ce, every plastic
     ForwardEuler sub-step has dsigma = Ce : de_e, so the whole sub-stepped
     increment must satisfy dsigma = Ce : de_e to round-off. Ce is read as the
     TanType-1 tangent of a zero-increment hold step just before. Holds
     whether or not ForwardEuler's multiplier is right (WP-158), and breaks on
     either defect: garbage moduli (dsigma != Ce : de_e) or an elastic strain
     that is never advanced (de_e = one sub-step's).
  2. ``test_energy_halves_equal_the_committed_split`` [0, 6] -- MaxEnergyInc's
     two halves equal the same increment as two COMMITTED steps of the
     integrator they call, up to the moduli each commit refreshes. (A split
     comparison does NOT work for the ForwardEuler schemes on a build without
     WP-158: there each committed FE step lands inside the yield surface and
     the next one re-intersects it elastically; sub-steps never do.)
  3. ``test_sub_stepped_schemes_are_deterministic`` [0, 4, 6, 7, 8, 9] -- a
     GUARD, not the regression gate: the same increment, re-run after a
     different deck, is bit-identical. Pre-fix it PASSES on the dev box (the
     stack garbage read as nG, nK was 0 there -- which is also why the
     sub-stepped stress did not move at all); WP-129 measured run-to-run
     differences for IntScheme 4 on another binary.

MEASURED 2026-10-02, Windows, Ladruno_scripts\\build.bat. Pre-fix = origin/ladruno
d63f49750 (SRC/material identical to the 117f56060 base). Fixed = WP-160:

  gate 1  identity gap        4, 7, 8, 9   pre 1.0 (dsigma = 0)   fixed 2.1e-14, 1.0e-14 x3
  gate 2  vs committed split  0            pre 1.0 (dsigma = 0)   fixed 4.4e-4
                              6  (e_e)     pre 0.45               fixed 6.8e-15
  gate 3  re-run              all          pre True               fixed True
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a, pytest.mark.t0m]

# Toyoura reference set (the WP-110 tangent gate's MATP), medium-dense.
MATP = dict(G0=125.0, nu=0.05, e_init=0.75, Mc=1.25, c=0.712, lambda_c=0.019,
            e0=0.934, ksi=0.7, P_atm=101.3, m=0.01, h0=7.05, Ch=0.968, nb=1.1,
            A0=0.704, nd=3.5, z_max=4.0, cz=600.0, Rho=1.7)
_ORDER = ["G0", "nu", "e_init", "Mc", "c", "lambda_c", "e0", "ksi", "P_atm",
          "m", "h0", "Ch", "nb", "A0", "nd", "z_max", "cz", "Rho"]
COORDS = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
          5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}

IDENTITY_RTOL = 1e-10    # fixed 1-2e-14; pre-fix 1.0 (dsigma == 0)
SPLIT_RTOL = 5e-3        # fixed 4.4e-4 (0), 7e-15 (6); pre-fix 1.0 (0), 0.45 (6, e_e)


def _set_scheme(scheme, ptag=1):
    """Runtime switch; argv[1] of ManzariDafalias::setParameter is the material tag."""
    ops.parameter(ptag, "element", 1, "IntegrationScheme", 1)
    ops.updateParameter(ptag, float(scheme))
    ops.remove("parameter", ptag)


def _history():
    """Isotropic compression (elastic, 5 steps), then a deviatoric push with
    shear (40 steps); returns the strain points and the last step's increment
    (largest component 1e-4: MaxStrainInc takes 11 sub-steps of it)."""
    pts = [np.eye(3) * -1.5e-3 * (i + 1) / 5 for i in range(5)]
    dev = np.array([[0.5, 0.2, 0.0], [0.2, 0.5, 0.1], [0.0, 0.1, -1.0]]) * 4e-3
    for i in range(1, 41):
        pts.append(pts[4] + dev * i / 40)
    return pts, dev / 40


def run(strains, last_scheme, n_last=1):
    """One SSPbrick through 3x3 strains (zero free DOFs, so the Gauss point
    sees exactly these): IntScheme 1 throughout except the LAST `n_last`
    steps. Returns {stress, alpha, estrains} before and after those steps."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in COORDS.items():
        ops.node(tag, float(x), float(y), float(z))
    ops.fix(1, 1, 1, 1)
    p = MATP
    ops.nDMaterial("ManzariDafalias", 1, *[p[k] for k in _ORDER], 1, 1, 1, 1e-7, 1e-7)
    ops.element("SSPbrick", 1, *range(1, 9), 1)
    ops.constraints("Transformation")    # 'Plain' drops non-homogeneous SPs
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-9, 5, 0)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    times = [float(t) for t in range(len(strains) + 1)]
    pts = [np.zeros((3, 3))] + list(strains)
    tag = 100
    for n in range(2, 9):
        X = np.array(COORDS[n], dtype=float)
        for d in range(3):
            tag += 1
            ops.timeSeries("Path", tag, "-time", *times,
                           "-values", *[float((E @ X)[d]) for E in pts])
            ops.pattern("Plain", tag, tag)
            ops.sp(n, d + 1, 1.0)

    def snap():
        out = {k: np.array(ops.eleResponse(1, k)) for k in ("stress", "alpha", "estrains")}
        out["tangent"] = np.array(ops.eleResponse(1, "tangent")).reshape(6, 6)
        return out

    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for i in range(len(strains)):
        if i == 5:
            ops.updateMaterialStage("-material", 1, "-stage", 1)
        if i == len(strains) - n_last:
            before = snap()
            if last_scheme != 1:
                _set_scheme(last_scheme)
        assert ops.analyze(1) == 0, "prescribed step failed"
    return before, snap()


# sub-steps each scheme takes of the last increment (largest component 1e-4)
N_SUB = {7: 11, 8: 11, 9: 11, 4: 2, 0: 2, 6: 2}   # floor(1e-4/1e-5)+1 ; MaxEnergyInc halves
SINGLE = {7: 5, 8: 5, 9: 5, 4: 5, 0: 1, 6: 3}     # the integrator each sub-step calls


def increments(last_scheme, n_split=1, hold=False):
    """The last increment under `last_scheme`, taken in ONE step or split into
    `n_split` equal COMMITTED steps. `hold`: a zero-increment step first, so
    the "before" tangent is the committed Ce = GetStiffness(mK, mG) -- an
    elastic step's TanType-1 tangent -- which is what a sub-stepper hands its
    sub-steps after WP-160."""
    pts, step = _history()
    head = pts[:-1] + ([pts[-2]] if hold else [])
    tail = [head[-1] + step * (j + 1) / n_split for j in range(n_split)]
    b, a = run(head + tail, last_scheme, n_last=n_split)
    if hold:
        return {k: a[k] - b[k] for k in ("stress", "alpha", "estrains")}, a, b["tangent"]
    return {k: a[k] - b[k] for k in ("stress", "alpha", "estrains")}, a


def rel(x, ref):
    return float(np.linalg.norm(x - ref) / np.linalg.norm(ref))


def elastic_identity_gap(scheme):
    """|dsigma - Ce : de_e| / max(|dsigma|, |Ce : de_e|) over the sub-stepped
    increment. With the moduli frozen at the committed Ce, every plastic
    ForwardEuler sub-step has dsigma = Ce : (de - dgamma R) and de_e = de -
    dgamma R, so the SUM holds exactly -- whether or not ForwardEuler's
    multiplier is right (WP-158). Garbage moduli break it (dsigma != Ce:de_e),
    and so does an elastic strain that is never advanced (de_e = one sub-step's)."""
    d, a, Ce = increments(scheme, hold=True)
    ce = Ce @ d["estrains"]
    return float(np.linalg.norm(d["stress"] - ce) / max(np.linalg.norm(d["stress"]), np.linalg.norm(ce))), d


# ------------------------------------------------------------------ gates ---
@pytest.mark.parametrize("scheme", [4, 7, 8, 9])
def test_forward_euler_substeps_keep_the_elastic_bookkeeping(scheme):
    gap, d = elastic_identity_gap(scheme)
    assert np.linalg.norm(d["stress"]) > 0.0, f"IntScheme {scheme}: the increment did not move the stress"
    assert gap < IDENTITY_RTOL, (
        f"IntScheme {scheme}: dsigma != Ce_committed : de_e over the sub-stepped increment "
        f"(relative gap {gap:.3e})", d)


@pytest.mark.parametrize("scheme", [0, 6])
def test_energy_halves_equal_the_committed_split(scheme):
    """MaxEnergyInc's two halves vs the same increment as two COMMITTED steps of
    the integrator it calls (ModifiedEuler for 0, RungeKutta4 for 6); they
    differ only by the moduli the commit refreshes (RK4 re-evaluates them
    itself, so 6 agrees to round-off)."""
    d, a = increments(scheme)
    dr, ar = increments(SINGLE[scheme], N_SUB[scheme])
    d1, a1 = increments(SINGLE[scheme])
    assert any(not np.array_equal(a[k], a1[k]) for k in d), (
        f"IntScheme {scheme} did not sub-step -- the gate would be vacuous")
    for k in ("stress", "alpha", "estrains"):
        r = rel(d[k], dr[k])
        assert r < SPLIT_RTOL, (
            f"IntScheme {scheme}: the halved {k} increment differs from two committed "
            f"IntScheme-{SINGLE[scheme]} steps by {r:.3e}", d[k], dr[k])


def _dirty_the_stack():
    """A different deck in between: other values left where nG, nK used to be read."""
    pts, step = _history()
    run([E * 3.0 for E in pts], 2)


@pytest.mark.parametrize("scheme", [0, 4, 6, 7, 8, 9])
def test_sub_stepped_schemes_are_deterministic(scheme):
    _, first = increments(scheme)
    _dirty_the_stack()
    _, again = increments(scheme)
    for k in ("stress", "alpha", "estrains"):
        assert np.array_equal(first[k], again[k]), (
            f"IntScheme {scheme}: the same increment gave a different {k} on a re-run",
            first[k], again[k])


if __name__ == "__main__":
    for s_ in (4, 7, 8, 9, 0, 6):
        gap, d = elastic_identity_gap(s_)
        print("IntScheme %d: elastic identity gap %.3e  |dsigma| %.4e" % (s_, gap, np.linalg.norm(d["stress"])))
    for s_ in (0, 6, 4, 9):
        d, a = increments(s_)
        dr, _ = increments(SINGLE[s_], N_SUB[s_])
        print("IntScheme %d vs %d committed IntScheme-%d steps: stress %.3e  alpha %.3e  estrains %.3e" % (
            s_, N_SUB[s_], SINGLE[s_], rel(d["stress"], dr["stress"]), rel(d["alpha"], dr["alpha"]),
            rel(d["estrains"], dr["estrains"])))
    for s_ in (0, 4, 6, 7, 8, 9):
        _, a = increments(s_)
        _dirty_the_stack()
        _, again = increments(s_)
        print("IntScheme %d re-run bit-identical: %s" % (s_, all(np.array_equal(a[k], again[k]) for k in ("stress", "alpha", "estrains"))))
