"""EnergyBalanceRecorder (per-region energy sidecar, RECORDER classTag 33000) —
Zone-A T1 check.

Zone-A port of the EnergyBalance check (test 9) from
``Ladruno_scripts/_verify_explicit.py`` (band-remediated 26 -> 33000). The
recorder writes a per-step energy ledger; with ``-time`` the columns are

    time  KE  IE  DW  ULW  RES  ERR%

(confirmed against SRC/recorder/EnergyBalanceRecorder.cpp, column header
{"KE","IE","DW","ULW","RES","ERR%"}, time prepended when ``-time`` is given).

The primary check is the getMass()/getDamp() shared-``theMatrix`` aliasing-fix
regression: with ELEMENT mass (Truss ``-rho``) the recorder must read kinetic
energy from the element MASS matrix, not the (zeroed) damping matrix. Before the
fix the element KE came out ~0; afterwards it tracks 0.5 * M * v^2. This recorder
already runs green on Zone-A CI (CDL-10 exercises the same KE/IE columns).
"""
import math
import re

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _v2_cols(capfd):
    """Column names from the recorder's echoed header (the documented way to
    parse -v2 output — 'parse by column names, never by position'). Reads the
    output captured since the last readouterr(), so call once right after the
    run whose layout you need."""
    cap = capfd.readouterr()
    text = cap.out + cap.err
    m = re.findall(r"columns =( time)? \[model: ([^\]]+)\]", text)
    assert m, "no EnergyBalanceRecorder column echo captured"
    has_time, names = m[-1]
    return (["time"] if has_time else []) + names.split()


def build_bar(N, L, E, A, rho):
    """Fixed-free chain of N truss elements along x; returns element length."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    le = L / N
    for i in range(N + 1):
        ops.node(i + 1, i * le, 0.0)
    ops.fix(1, 1, 1)
    for i in range(1, N + 1):
        ops.fix(i + 1, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    for i in range(N):
        ops.element("Truss", i + 1, i + 1, i + 2, A, 1, "-rho", rho)
    return le


# ==========================================================================
# element-mass KE: the getMass()/getDamp() aliasing-fix regression (test 9)
# ==========================================================================
def test_energybalance_element_mass_ke(tmp_path):
    N = 20
    L = 10.0
    E = 100.0
    A = 1.0
    rho = 1.0
    v0 = 2.0
    efile = str(tmp_path / "bar_energy.txt")
    build_bar(N, L, E, A, rho)
    for i in range(1, N + 1):                 # initial x-velocity on every free node
        ops.setNodeVel(i + 1, 1, v0, "-commit")
    ops.recorder("EnergyBalance", "-file", efile, "-time")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("ExplicitBathe", 0.54)
    ops.analysis("Transient")
    ops.analyze(2, 1e-4)
    ops.wipe()                                # flush/close the recorder

    d = np.loadtxt(efile)
    d = np.atleast_2d(d)
    ke0 = d[0, 1]                             # col 0 = time, col 1 = KE
    M_total = rho * A * L                     # total bar mass (element-based)
    ke_expected = 0.5 * M_total * v0 ** 2     # ~ (fixed end excludes ~rho*A*le/2)
    assert ke0 > 0.5 * ke_expected, (
        "element-mass KE=%.3f vs ~0.5*M*v^2=%.3f — before the aliasing fix this "
        "was ~0 (read from the zeroed damping matrix)" % (ke0, ke_expected)
    )


# ==========================================================================
# column layout: with -time, exactly 7 columns (time + KE IE DW ULW RES ERR%)
# ==========================================================================
def _run_lnvd_bar(tmp_path, fname, v2, alpha):
    """Free-vibration bar with ExplicitBathe [-lnvd alpha]; returns the data."""
    N = 10
    L = 1.0
    E = 100.0
    A = 1.0
    rho = 1.0
    efile = str(tmp_path / fname)
    le = build_bar(N, L, E, A, rho)
    for i in range(1, N + 1):
        ops.setNodeVel(i + 1, 1, 1.0, "-commit")
    # NOTE (ADR-69): integrator BEFORE recorder — the -v2 channel columns are
    # fixed at the recorder's first record from the registry's declared set.
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    if alpha > 0.0:
        ops.integrator("ExplicitBathe", 0.54, "-lnvd", alpha)
    else:
        ops.integrator("ExplicitBathe", 0.54)
    if v2:
        ops.recorder("EnergyBalance", "-file", efile, "-time", "-v2")
    else:
        ops.recorder("EnergyBalance", "-file", efile, "-time")
    ops.analysis("Transient")
    c = math.sqrt(E / rho)
    ops.analyze(400, 0.2 * le / c)
    ops.wipe()
    d = np.loadtxt(efile)
    return np.atleast_2d(d)


# ==========================================================================
# ADR-69 v2 column geometry, order-INDEPENDENT (channel declarations are
# process-sticky, so earlier test FILES in the battery — e.g. the adr52
# ExplicitBathe-collapse tests constructing -lnvd — may have declared
# channels before this file runs; those then appear as benign zero columns).
# Fixed geometry: time | KE_ele KE_nod IE DW_ele DW_nod ULW | [channels...]
# | RES ERR%  -> cols 1..6 fixed, RES/ERR always the LAST TWO, channels (if
# any) sit at 7..ncol-3 with E_lnvd LAST among them (E_inject, E_lnvd order).
# ==========================================================================
def test_energybalance_v2_layout_no_producers(tmp_path):
    d = _run_lnvd_bar(tmp_path, "v2_plain.txt", v2=True, alpha=0.0)
    ncol = d.shape[1]
    assert ncol >= 9, (
        "expected >=9 columns (time + 6 split + RES ERR%%); got %d" % ncol
    )
    # THIS run has no producer: any channel columns (declared by earlier test
    # files in the same process) must be identically zero.
    if ncol > 9:
        assert np.allclose(d[:, 7:ncol - 2], 0.0, atol=1e-14), (
            "no-producer run: channel columns must be all-zero"
        )
    ke_ele, ke_nod = d[0, 1], d[0, 2]
    assert ke_ele > 0.0 and abs(ke_nod) < 1e-12, (
        "element-mass bar: KE must sit in KE_ele (got ele=%g nod=%g)" % (ke_ele, ke_nod)
    )


# ==========================================================================
# ADR-69 v2: E_lnvd channel closes what v1 leaks into RES. Free vibration +
# FLAC local damping: legacy RES drifts by exactly the (unaccounted) LNVD
# dissipation; v2 publishes it (E_lnvd column) and RES stays put.
# ==========================================================================
def test_energybalance_v2_lnvd_closure(tmp_path, capfd):
    alpha = 0.4
    d_legacy = _run_lnvd_bar(tmp_path, "lnvd_legacy.txt", v2=False, alpha=alpha)
    capfd.readouterr()   # drop the legacy run's echo
    d_v2 = _run_lnvd_bar(tmp_path, "lnvd_v2.txt", v2=True, alpha=alpha)
    cols = _v2_cols(capfd)

    # legacy: time KE IE DW ULW RES ERR% -> RES col 5
    res_l = d_legacy[:, 5]
    drift_legacy = res_l[-1] - res_l[0]

    # v2 with LNVD declared: locate E_lnvd BY NAME from the echoed header.
    # (Back-indexing d[:, -3] broke the moment a later channel — E_modal,
    # ADR-69 P2 — could be declared by an earlier test in the same process;
    # name lookup is the documented parse contract and is order-proof.)
    assert d_v2.shape[1] >= 10, (
        "expected >=10 columns (time + 6 split + E_lnvd + RES ERR%%); got %d"
        % d_v2.shape[1]
    )
    e_lnvd = d_v2[:, cols.index("E_lnvd")]
    res_v2 = d_v2[:, cols.index("RES")]
    drift_v2 = res_v2[-1] - res_v2[0]

    assert e_lnvd[-1] > 0.0, "LNVD dissipation must accumulate positive work"
    assert np.all(np.diff(e_lnvd) > -1e-14), "E_lnvd must be monotone"
    # the legacy RES drift IS the unaccounted LNVD dissipation
    assert drift_legacy > 0.0
    assert 0.7 * drift_legacy < e_lnvd[-1] < 1.3 * drift_legacy, (
        "published E_lnvd=%g should match the legacy RES drift=%g"
        % (e_lnvd[-1], drift_legacy)
    )
    # closure: v2 RES drift shrinks by at least 5x
    assert abs(drift_v2) < 0.2 * abs(drift_legacy), (
        "v2 RES drift %g should be <<%g (legacy)" % (drift_v2, drift_legacy)
    )


# ==========================================================================
# ADR-69 P1.5: ASDAbsorbingBoundary2D bottom compliant-base injection closes
# the same way LysmerTriangle's does. addBaseActions() is a pure external
# source (time-series velocity * fixed coefficients on the SOIL-side nodes,
# mutually exclusive with the lateral free-field mechanism via BND_BOTTOM)
# that leaks into the legacy IE integral; -v2 rebuckets it into E_inject.
# Ghost nodes are explicitly fixed (belt-and-suspenders on top of the
# element's own internal penalty pin) to keep the system well-posed.
# ==========================================================================
def _build_asd_column(amp_scale, n_layers=3):
    E, nu, rho = 2.0 * 2000.0 * 200.0 ** 2 * 1.25, 0.25, 2000.0
    G = 2000.0 * 200.0 ** 2
    dt = 1.0e-3
    nstep = 60
    t = np.arange(nstep + 1) * dt
    f0, t0, amp = 10.0, 0.02, 0.05
    a = (np.pi * f0 * (t - t0)) ** 2
    vals = (amp * amp_scale * (1.0 - 2.0 * a) * np.exp(-a)).tolist()

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1001, 0.0, -1.0)
    ops.node(1002, 1.0, -1.0)
    ops.fix(1001, 1, 1)
    ops.fix(1002, 1, 1)
    for row in range(n_layers + 1):
        ops.node(1 + 2 * row, 0.0, float(row))
        ops.node(2 + 2 * row, 1.0, float(row))
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    for row in range(n_layers):
        n1, n2 = 1 + 2 * row, 2 + 2 * row
        n3, n4 = 2 + 2 * (row + 1), 1 + 2 * (row + 1)
        ops.element("quad", 200 + row, n1, n2, n3, n4, 1.0, "PlaneStrain", 1)
    ops.timeSeries("Path", 1, "-dt", dt, "-values", *vals)
    ops.element("ASDAbsorbingBoundary2D", 101, 1001, 1002, 2, 1,
               G, nu, rho, 1.0, "B", "-fx", 1)
    ops.setParameter("-val", 1, "-ele", 101, "stage")
    return dt, nstep


def _run_asd_column(tmp_path, fname, v2, amp_scale):
    dt, nstep = _build_asd_column(amp_scale)
    efile = str(tmp_path / fname)
    rec = ["EnergyBalance", "-file", efile, "-time"]
    if v2:
        rec.append("-v2")
    ops.recorder(*rec)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 10)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    for _ in range(nstep):
        assert ops.analyze(1, dt) == 0
    ops.wipe()
    return np.atleast_2d(np.loadtxt(efile))


def test_energybalance_v2_asd_absorbing_injection_closure(tmp_path):
    d_legacy = _run_asd_column(tmp_path, "asd_legacy.txt", v2=False, amp_scale=1.0)
    d_v2 = _run_asd_column(tmp_path, "asd_v2.txt", v2=True, amp_scale=1.0)
    d_ctl = _run_asd_column(tmp_path, "asd_ctl.txt", v2=True, amp_scale=0.0)

    ie_legacy = d_legacy[-1, 2]
    eref_legacy = max(np.abs(d_legacy[:, 1:5]).max(), 1e-30)
    assert abs(ie_legacy) > 0.05 * eref_legacy, (
        "expected the legacy IE leak to be a significant fraction of E_ref; "
        "got IE=%g E_ref=%g" % (ie_legacy, eref_legacy)
    )

    # This model's ASDAbsorbingBoundary2D unconditionally declares ABSORB_LEAK
    # in setDomain, so chInject is guaranteed true and E_inject sits at the
    # FIXED front-anchored offset col 7 (right after ULW at col 6) -
    # regardless of whether an earlier test in this process (or battery run)
    # also declared LNVD_WORK, which the recorder always writes AFTER
    # E_inject (see EnergyBalanceRecorder.cpp record(): chInject column is
    # written before chLnvd). Do NOT locate E_inject from the back like the
    # LNVD test does - E_lnvd (not E_inject) is what's guaranteed last.
    assert d_v2.shape[1] >= 10
    ie_v2 = d_v2[-1, 3]
    e_inject = d_v2[-1, 7]
    eref_v2 = max(np.abs(d_v2[:, 1:8]).max(), 1e-30)
    assert abs(ie_v2) < 0.05 * eref_v2, (
        "v2 IE should be near-truthful once the leak is rebucketed; got %g "
        "(E_ref %g)" % (ie_v2, eref_v2)
    )

    diff = ie_v2 - ie_legacy
    assert abs(e_inject - diff) < 0.1 * max(abs(diff), 1e-30), (
        "E_inject (%g) should match IE_v2 - IE_legacy (%g) - the same leak "
        "measured through two independent paths" % (e_inject, diff)
    )

    e_inject_ctl = d_ctl[-1, 7]
    assert abs(e_inject_ctl) < 1e-9, (
        "zero-amplitude control: E_inject should be ~0, got %g" % e_inject_ctl
    )


# ==========================================================================
# ADR-69 P1.6: ASDAbsorbingBoundary2D LATERAL free-field transfer closes the
# same way. addRffToSoil() (mutually exclusive with the bottom mechanism via
# BND_BOTTOM) writes the free-field column's transfer force onto the SOIL
# nodes only; it leaks into the legacy IE integral and -v2 rebuckets it into
# E_inject. Model: soil quad column at x in [0,1] + free-field node column
# at x=-1 with one "L" boundary element per layer. Upper FF nodes MUST stay
# free (stage-1 penalty skips vertical boundaries; the FF dynamics live in
# addKff/addMff/addCff). "L" takes no -fx/-fy (parser rejects them for
# non-bottom).
#
# Drive is INITIAL VELOCITY on both DOFs of every free node (cross-axis
# transfer pairing: mu-term soil-y force from FF du_x, lambda-term soil-x
# force from FF du_y — both axes must move for first-order work). NOT
# UniformExcitation, for two source-verified reasons: (1) the ASD element's
# addInertiaLoadToUnbalance is a deliberate no-op, so the FF masses are
# never driven and the FF column stays rigid; (2) the quads store -M*ug''
# in their element load vector Q and getResistingForce returns K*u-Q, so
# seismic input work hides in IE (IE = -DW exactly, RES accidentally
# closed) and would drown the lateral-transfer signal the -v2 rebucket
# targets. Mass-proportional Rayleigh (alphaM only) makes the FF column
# decay too — the lateral addClk dashpot writes only SOIL rows (one-way
# coupling), so an undamped FF column rings forever and its genuine spring
# energy would keep IE_end from settling.
# ==========================================================================
def _build_asd_lateral_column(amp_scale, n_layers=3):
    E, nu, rho = 2.0 * 2000.0 * 200.0 ** 2 * 1.25, 0.25, 2000.0
    G = 2000.0 * 200.0 ** 2
    dt = 1.0e-3
    nstep = 300
    v0 = 0.5 * amp_scale

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for row in range(n_layers + 1):
        ops.node(1 + 2 * row, 0.0, float(row))
        ops.node(2 + 2 * row, 1.0, float(row))
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    # free-field column at x=-1; ONLY the base fixed
    for row in range(n_layers + 1):
        ops.node(1000 + row, -1.0, float(row))
    ops.fix(1000, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    for row in range(n_layers):
        n1, n2 = 1 + 2 * row, 2 + 2 * row
        n3, n4 = 2 + 2 * (row + 1), 1 + 2 * (row + 1)
        ops.element("quad", 200 + row, n1, n2, n3, n4, 1.0, "PlaneStrain", 1)
    for row in range(n_layers):
        ops.element("ASDAbsorbingBoundary2D", 100 + row,
                    1000 + row, 1000 + row + 1,
                    1 + 2 * row, 1 + 2 * (row + 1),
                    G, nu, rho, 1.0, "L")
    for row in range(n_layers):
        ops.setParameter("-val", 1, "-ele", 100 + row, "stage")
    for row in range(1, n_layers + 1):
        ops.setNodeVel(1 + 2 * row, 1, v0, "-commit")
        ops.setNodeVel(1 + 2 * row, 2, v0, "-commit")
        ops.setNodeVel(2 + 2 * row, 1, v0, "-commit")
        ops.setNodeVel(2 + 2 * row, 2, v0, "-commit")
        ops.setNodeVel(1000 + row, 1, v0, "-commit")
        ops.setNodeVel(1000 + row, 2, v0, "-commit")
    ops.rayleigh(20.0, 0.0, 0.0, 0.0)
    return dt, nstep


def _run_asd_lateral_column(tmp_path, fname, v2, amp_scale):
    dt, nstep = _build_asd_lateral_column(amp_scale)
    efile = str(tmp_path / fname)
    rec = ["EnergyBalance", "-file", efile, "-time"]
    if v2:
        rec.append("-v2")
    ops.recorder(*rec)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 10)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    for _ in range(nstep):
        assert ops.analyze(1, dt) == 0
    ops.wipe()
    return np.atleast_2d(np.loadtxt(efile))


def test_energybalance_v2_asd_lateral_freefield_closure(tmp_path):
    d_legacy = _run_asd_lateral_column(tmp_path, "asdlat_legacy.txt", v2=False,
                                       amp_scale=1.0)
    d_v2 = _run_asd_lateral_column(tmp_path, "asdlat_v2.txt", v2=True,
                                   amp_scale=1.0)
    d_ctl = _run_asd_lateral_column(tmp_path, "asdlat_ctl.txt", v2=True,
                                    amp_scale=0.0)

    # E_inject at the FIXED front-anchored offset col 7 (right after ULW) -
    # this model's boundary elements unconditionally declare ABSORB_LEAK, and
    # the recorder writes E_inject before E_lnvd regardless of declaration
    # order elsewhere in the process. Never index channels from the back.
    assert d_v2.shape[1] >= 10
    ie_legacy = d_legacy[-1, 2]
    ie_v2 = d_v2[-1, 3]
    e_inject = d_v2[-1, 7]

    # the closure identity is exact by kernel algebra (IE_display = IE_raw - L,
    # E_inject = -L): the same leak measured through two independent paths
    diff = ie_v2 - ie_legacy
    assert abs(e_inject) > 0.0, "lateral transfer published no energy"
    assert abs(e_inject - diff) < 0.1 * max(abs(diff), 1e-30), (
        "E_inject (%g) should match IE_v2 - IE_legacy (%g)" % (e_inject, diff)
    )

    # the leak is a first-order fraction of the run's energy scale, and the
    # rebucketed IE is truthful (small relative to the same scale)
    eref = max(np.abs(d_v2[:, 1:8]).max(), 1e-30)
    assert abs(ie_legacy) > 0.05 * eref, (
        "expected a significant legacy lateral-transfer leak; got IE=%g "
        "(E_ref %g)" % (ie_legacy, eref)
    )
    assert abs(ie_v2) < 0.05 * eref, (
        "v2 IE should be near-truthful after the rebucket; got %g (E_ref %g)"
        % (ie_v2, eref)
    )

    e_inject_ctl = d_ctl[-1, 7]
    assert abs(e_inject_ctl) < 1e-9, (
        "zero-amplitude control: E_inject should be ~0, got %g" % e_inject_ctl
    )


# ==========================================================================
# ADR-69 P2: E_modal channel closes what v1 leaks into RES. Modal damping is
# applied by the integrator straight into the SOE B vector (-C_modal*v via
# setB) - it never touches an element/node damping matrix, so the recorder's
# DW path is structurally blind to it: legacy RES drifts by exactly the
# unaccounted modal dissipation. IncrementalIntegrator::commit() publishes
# the trapezoid of v'C_modal v (captured at the last converged formUnbalance
# iterate) to MODAL_WORK; -v2 shows it as E_modal and RES stays put. Use
# modalDamping (matrix-consistent), NOT modalDampingQ - Q applies the force
# off the lagged velocity with no tangent term and can ANTI-damp light free
# vibration. FullGeneral: the modal damping matrix is dense.
# ==========================================================================
def _run_modal_bar(tmp_path, fname, v2, zeta):
    N = 20
    L = 10.0
    E = 100.0
    A = 1.0
    rho = 1.0
    v0 = 2.0
    efile = str(tmp_path / fname)
    build_bar(N, L, E, A, rho)
    for i in range(1, N + 1):
        ops.setNodeVel(i + 1, 1, v0, "-commit")
    ops.eigen(1)
    ops.modalDamping(zeta)
    rec = ["EnergyBalance", "-file", efile, "-time"]
    if v2:
        rec.append("-v2")
    ops.recorder(*rec)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-10, 10)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    c = math.sqrt(E / rho)
    le = L / N
    assert ops.analyze(400, 0.2 * le / c) == 0
    ops.wipe()
    return np.atleast_2d(np.loadtxt(efile))


def test_energybalance_v2_modal_closure(tmp_path, capfd):
    zeta = 0.05
    d_legacy = _run_modal_bar(tmp_path, "modal_legacy.txt", v2=False, zeta=zeta)
    capfd.readouterr()
    d_v2 = _run_modal_bar(tmp_path, "modal_v2.txt", v2=True, zeta=zeta)
    cols = _v2_cols(capfd)

    # legacy: time KE IE DW ULW RES ERR% -> RES col 5 drifts by the
    # unaccounted modal dissipation
    res_l = d_legacy[:, 5]
    drift_legacy = res_l[-1] - res_l[0]
    assert drift_legacy > 0.0, (
        "legacy RES should drift positive by the modal dissipation; got %g"
        % drift_legacy
    )

    assert "E_modal" in cols, "E_modal column missing from -v2 layout"
    e_modal = d_v2[:, cols.index("E_modal")]
    res_v2 = d_v2[:, cols.index("RES")]
    drift_v2 = res_v2[-1] - res_v2[0]

    assert e_modal[-1] > 0.0, "modal dissipation must accumulate positive work"
    assert np.all(np.diff(e_modal) > -1e-9), (
        "E_modal must be monotone (trapezoid of v'Cv >= 0 rates)"
    )
    # the legacy RES drift IS the unaccounted modal dissipation
    assert 0.7 * drift_legacy < e_modal[-1] < 1.3 * drift_legacy, (
        "published E_modal=%g should match the legacy RES drift=%g"
        % (e_modal[-1], drift_legacy)
    )
    # closure: v2 RES drift shrinks by at least 5x
    assert abs(drift_v2) < 0.2 * abs(drift_legacy), (
        "v2 RES drift %g should be <<%g (legacy)" % (drift_v2, drift_legacy)
    )


# ==========================================================================
# ADR-69 P2: E_hg hourglass attribution (mechanism C recorder-pull). The URI
# LadrunoBrick's stabilization force is inside getResistingForce, so its
# work is booked in IE as if it were physical strain energy; the recorder
# probes elements for the "hourglassEnergy" response at initialize and
# rebuckets IE_display = IE_raw - E_hg with E_hg its own column. E_hg is the
# INSTANTANEOUS stored stabilization energy for the stiffness flavour (can
# oscillate - no monotonicity assert), and the rebucket cancels in RES by
# kernel algebra, so E_hg == IE_legacy - IE_v2 exactly (deterministic
# reruns) and RES is invariant. Dynamic bending cantilever from an initial
# transverse tip velocity excites the hourglass modes.
# ==========================================================================
def _run_hg_cantilever(tmp_path, fname, v2, form_args, nx=4):
    L, E, nu, rho = 10.0, 1000.0, 0.0, 1.0
    efile = str(tmp_path / fname)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    dx = L / nx

    def nt(i, j, k):
        return 1 + i + (nx + 1) * j + (nx + 1) * 2 * k

    for k in range(2):
        for j in range(2):
            for i in range(nx + 1):
                ops.node(nt(i, j, k), i * dx, j * 1.0, k * 1.0)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    for i in range(nx):
        ops.element("LadrunoBrick", 1 + i,
                    nt(i, 0, 0), nt(i + 1, 0, 0), nt(i + 1, 1, 0), nt(i, 1, 0),
                    nt(i, 0, 1), nt(i + 1, 0, 1), nt(i + 1, 1, 1), nt(i, 1, 1),
                    1, *form_args)
    for k in range(2):
        for j in range(2):
            ops.fix(nt(0, j, k), 1, 1, 1)
    # initial transverse velocity, linear in x (bending-dominated start)
    for k in range(2):
        for j in range(2):
            for i in range(1, nx + 1):
                ops.setNodeVel(nt(i, j, k), 3, 0.5 * i / nx, "-commit")
    # precision 12: the attribution identity below compares E_hg against the
    # NEAR-CANCELLING difference IE_legacy - IE_v2; at the default 6 sig figs
    # the file quantization (~1e-7 of IE) dwarfs the identity tolerance.
    rec = ["EnergyBalance", "-file", efile, "-time", "-precision", 12]
    if v2:
        rec.append("-v2")
    ops.recorder(*rec)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-10, 10)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    assert ops.analyze(150, 2.0e-3) == 0
    ops.wipe()
    return np.atleast_2d(np.loadtxt(efile))


def test_energybalance_v2_hourglass_attribution(tmp_path, capfd):
    uri = ["-formulation", "uri", "-hourglass", "stiffness"]
    d_legacy = _run_hg_cantilever(tmp_path, "hg_legacy.txt", v2=False,
                                  form_args=uri)
    capfd.readouterr()
    d_v2 = _run_hg_cantilever(tmp_path, "hg_v2.txt", v2=True, form_args=uri)
    cols = _v2_cols(capfd)

    assert "E_hg" in cols, "E_hg column missing from -v2 layout"
    e_hg = d_v2[:, cols.index("E_hg")]
    ie_v2 = d_v2[:, cols.index("IE")]
    res_v2 = d_v2[:, cols.index("RES")]

    # the URI cantilever genuinely hourglasses: nonzero stored stabilization
    # energy along the run (instantaneous - check the running max, not the
    # final sample, which may sit near a zero crossing)
    assert np.abs(e_hg).max() > 0.0, "URI bending run published no E_hg"

    # attribution identity: E_hg == IE_legacy - IE_v2 at every sample
    # (deterministic reruns of the same model; exact by kernel algebra)
    ie_legacy = d_legacy[:, 2]
    scale = max(np.abs(ie_legacy).max(), 1e-30)
    assert np.allclose(e_hg, ie_legacy - ie_v2, atol=1e-9 * scale), (
        "E_hg must equal IE_legacy - IE_v2 (the rebucket, measured through "
        "two independent paths)"
    )

    # RES invariance: the rebucket cancels in RES (attribution, not closure)
    res_legacy = d_legacy[:, 5]
    assert np.allclose(res_v2, res_legacy, atol=1e-9 * max(scale, 1.0)), (
        "rebucketing E_hg must not change RES"
    )

    # control: fully-integrated std twin -> the response exists (column
    # present) but the stored stabilization energy is identically zero
    capfd.readouterr()
    d_std = _run_hg_cantilever(tmp_path, "hg_std.txt", v2=True,
                               form_args=["-formulation", "std"])
    cols_std = _v2_cols(capfd)
    assert "E_hg" in cols_std
    e_hg_std = d_std[:, cols_std.index("E_hg")]
    assert np.abs(e_hg_std).max() < 1e-12, (
        "std formulation must publish zero hourglass energy"
    )


# ==========================================================================
# column layout: with -time, exactly 7 columns (time + KE IE DW ULW RES ERR%)
# ==========================================================================
def test_energybalance_time_column_layout(tmp_path):
    N = 5
    L = 1.0
    E = 100.0
    A = 1.0
    rho = 1.0
    efile = str(tmp_path / "cols.txt")
    le = build_bar(N, L, E, A, rho)
    ops.setNodeVel(N + 1, 1, 1.0, "-commit")
    ops.recorder("EnergyBalance", "-file", efile, "-time")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("ExplicitBathe", 0.54)
    ops.analysis("Transient")
    c = math.sqrt(E / rho)
    ops.analyze(5, 0.5 * le / c)
    ops.wipe()

    rows = []
    with open(efile) as fh:
        for line in fh:
            vals = line.split()
            if vals:
                rows.append(vals)
    assert rows, "EnergyBalance recorder produced no output"
    # time + model block (KE IE DW ULW RES ERR%) = 1 + 6 = 7 columns
    assert all(len(r) == 7 for r in rows), (
        "expected 7 columns per row (time + KE IE DW ULW RES ERR%%); got widths %s"
        % sorted({len(r) for r in rows})
    )
