"""Mass-scaling MOTIVATING-CASE validation (Ladruno ADR 37, Tier 2 / T-ZONEB) — Zone-B.

The Zone-A fidelity tests (test_massScaling_validation.py) prove SMS on tiny hand-built
chains. This is the honest "does it deliver its motivating value" proof on a *real*
gmsh-meshed 3D hex bar with a localized REFINEMENT ZONE of small elements — the actual
SSI / contact / pile use case, where a few tiny elements in an otherwise coarse model
choke the explicit step.

Bar: 1x1x10 prism (E=1000, nu=0, rho=1 -> clean axial speed c=sqrt(E/rho)), clamped at
z=0, meshed 2x2 across with a graded axial extrude: a short FINE band near the fixed end
(dz=0.25, 2x finer than the coarse dz=0.5 the cross-section already caps) and coarse
elsewhere. The fine band's per-element step (~0.25/c) is half the bulk step (~0.5/c), so
plain central difference at the bulk step is UNSTABLE — exactly what SMS is for.

The headline result this pins:
  * NECESSITY  — plain CentralDifferenceLadruno at the bulk dt DIVERGES (the tiny elements
    are unstable) while CentralDifferenceSMS at the same dt stays bounded. A no-op SMS
    (no injection) diverges here too -> caught.
  * FIDELITY   — placed in a LOW-participation zone (near the clamp), the fictitious mass
    (~20% of the model) barely moves the structural dynamics: the fundamental frequency
    shifts < 1% (eigen), and the free-end axial vibration period + peak amplitude track the
    fine-dt reference. (Period/peak isolate the mass-scaling fidelity; a raw full-history
    RMS would also fold in the bulk-dt time-integration dispersion, which is the *point* of
    running at the big step, not a scaling error.)
  * REPORT     — the achieved %added-mass and df1 are surfaced in the assert messages.

Zone-B: needs gmsh for meshing (raw gmsh API, not apeGmsh); runs nightly.
"""
import math

import pytest

from _testbed import ops

gmsh = pytest.importorskip("gmsh")
pytestmark = [pytest.mark.zone_b]

E, NU, RHO = 1000.0, 0.0, 1.0
C = math.sqrt(E / RHO)          # 1-D axial wave speed (nu=0)
A, L = 1.0, 10.0               # cross-section side, bar length (z axis)
COARSE_STEP = 0.5 / C          # bulk element-pencil step (cross-section 0.5 caps it)
DT_TARGET = 0.9 * COARSE_STEP  # SMS target: the bulk step (just under it)
SEED = 0.02                    # axial pull amplitude (linear ramp -> excites mode 1)
N_FINE = 16                    # 4 axial fine layers x 4 (2x2) cross-section hexes


def _mesh():
    """2x2-across prism, graded axial extrude: coarse(0->0.5) | FINE 4 layers(0.5->1.5,
    dz=0.25) | coarse(1.5->10). gmsh hex node order == OpenSees Brick order."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("bar")
        occ = gmsh.model.occ
        occ.addRectangle(0, 0, 0, A, A)
        occ.synchronize()
        for c in gmsh.model.getEntities(1):
            gmsh.model.mesh.setTransfiniteCurve(c[1], 3)     # 2 divisions -> 2x2
        gmsh.model.mesh.setTransfiniteSurface(1)
        gmsh.model.mesh.setRecombine(2, 1)
        occ.synchronize()
        # numElements per band over cumulative normalized heights
        occ.extrude([(2, 1)], 0, 0, L, numElements=[1, 4, 17],
                    heights=[0.5 / L, 1.5 / L, 1.0], recombine=True)
        occ.synchronize()
        gmsh.model.mesh.generate(3)
        nt, nc, _ = gmsh.model.mesh.getNodes()
        nodes = {int(t): (nc[3 * i], nc[3 * i + 1], nc[3 * i + 2]) for i, t in enumerate(nt)}
        hexes = []
        for et, _, en in zip(*gmsh.model.mesh.getElements(3)):
            if et == 5:                                      # 8-node hex
                en = [int(x) for x in en]
                hexes += [en[8 * k:8 * k + 8] for k in range(len(en) // 8)]
        return nodes, hexes
    finally:
        gmsh.finalize()


# Build once per SESSION, lazily, and reuse the connectivity (meshing is the
# slow part). LAZY on purpose — this used to run at MODULE level, i.e. during
# pytest COLLECTION of the full battery: gmsh.initialize() REPLACES the
# process's native Win32 PATH (os.environ keeps the original snapshot, so
# python code never notices), after which every later test that spawns a
# child process breaks — bare-name CreateProcess lookups (WinError 2) and
# MinGW-exe DLL resolution (0xC0000135) both search the NATIVE path. The
# conftest _native_path_resync guard heals the env between tests, but gmsh
# must still never run at import time. See LEDGER_quirks.
_MESH_CACHE = None


def _mesh_once():
    global _MESH_CACHE
    if _MESH_CACHE is None:
        _MESH_CACHE = _mesh()
    return _MESH_CACHE


def _build(seed=False):
    _NODES, _HEXES = _mesh_once()
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    used = set(n for hx in _HEXES for n in hx)
    for t in used:
        ops.node(int(t), *_NODES[t])
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    for i, hx in enumerate(_HEXES, 1):
        ops.element("LadrunoBrick", i, *[int(n) for n in hx], 1,
                    "-formulation", "uri", "-hourglass", "stiffness", "-lumped")
    zmax = max(_NODES[t][2] for t in used)
    for t in used:
        if abs(_NODES[t][2]) < 1e-9:
            ops.fix(int(t), 1, 1, 1)                         # clamp z=0 face
    tip = sorted(t for t in used if abs(_NODES[t][2] - zmax) < 1e-9)[0]
    if seed:
        for t in used:                                       # axial (z) ramp -> mode 1
            ops.setNodeDisp(int(t), 3, SEED * _NODES[t][2] / L, "-commit")
    return used, tip


def _run(integrator_args, dt, T, cap=1.0e3):
    """Seeded free axial vibration; tip z-history. Stops early if |u| blows past cap."""
    used, tip = _build(seed=True)
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    t = [ops.getTime()]
    u = [ops.nodeDisp(tip, 3)]
    diverged = False
    for _ in range(int(round(T / dt))):
        if ops.analyze(1, dt) != 0:
            break
        ui = ops.nodeDisp(tip, 3)
        if not math.isfinite(ui) or abs(ui) > cap:
            diverged = True
            break
        t.append(ops.getTime())
        u.append(ui)
    return t, u, diverged


def _dt_cr():
    _build()
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno", "-cfl")
    ops.analysis("Transient")
    ops.analyze(1, 1.0e-9)
    return ops.criticalTimeStep()


def _f1(after_sms=False):
    """Fundamental frequency, optionally after one SMS domainChanged injection."""
    _build()
    if after_sms:
        ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 1)
        ops.algorithm("Linear")
        ops.integrator("CentralDifferenceSMS", DT_TARGET)
        ops.analysis("Transient")
        ops.analyze(1, 0.5 * DT_TARGET)          # commits the fictitious-mass injection
    ops.system("FullGeneral")
    lam = ops.eigen("-fullGenLapack", 1)
    return math.sqrt(lam[0]) / (2.0 * math.pi)


def _period(t, u):
    um = sum(u) / len(u)
    c = [x - um for x in u]
    cr = []
    for i in range(1, len(c)):
        if c[i - 1] <= 0.0 < c[i]:
            fr = -c[i - 1] / (c[i] - c[i - 1])
            cr.append(t[i - 1] + fr * (t[i] - t[i - 1]))
    if len(cr) < 2:
        return None
    return sum(cr[k + 1] - cr[k] for k in range(len(cr) - 1)) / (len(cr) - 1)


def _added_mass_pct(capsys=None):
    """Run one SMS step with -verbose, return (added_pct, n_scaled, n_elem) from the report."""
    import io
    import contextlib
    _build()
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    buf = io.StringIO()
    with contextlib.redirect_stderr(buf):
        ops.integrator("CentralDifferenceSMS", DT_TARGET, "-verbose")
        ops.analysis("Transient")
        ops.analyze(1, 0.5 * DT_TARGET)
    import re
    err = buf.getvalue()
    mf = re.search(r"added mass ([0-9.eE+\-]+)% of total element \(translational\) mass", err)
    ms = re.search(r"scaled (\d+)/(\d+)", err)
    return (float(mf.group(1)) if mf else float("nan"),
            int(ms.group(1)) if ms else -1,
            int(ms.group(2)) if ms else -1)


# ==========================================================================
# T-ZONEB: SMS on a real gmsh-meshed refined 3D bar — necessity + fidelity + report
# ==========================================================================
def test_zoneb_refined_bar_sms_fidelity():
    _NODES, _HEXES = _mesh_once()
    assert len(_HEXES) > 40, "mesh too small (%d hexes)" % len(_HEXES)
    dtcr = _dt_cr()
    # the fine band genuinely governs the global step, well below the bulk target
    assert 0.0 < dtcr < 0.7 * DT_TARGET, (
        "fine-zone dt_cr=%.4g must be well below the bulk dtTarget=%.4g (else there is no "
        "refinement to scale)" % (dtcr, DT_TARGET)
    )

    # --- scaling actually happened, and ONLY the fine band was scaled ---
    added_pct, n_scaled, n_elem = _added_mass_pct()
    assert n_scaled == N_FINE, (
        "SMS should scale exactly the %d fine-band hexes; scaled %d/%d" % (N_FINE, n_scaled, n_elem)
    )
    assert added_pct > 5.0, (
        "added mass %.2f%% too small — scaling not meaningfully exercised" % added_pct
    )

    # --- NECESSITY: plain CD at the bulk dt diverges (tiny elements unstable); SMS holds.
    #     A no-op / broken SMS that injects nothing diverges here too -> this leg catches it.
    dt_bulk = 0.9 * DT_TARGET
    T = 3.0 * (4.0 * L / C)                       # ~3 fundamental periods
    _, u_cd, cd_div = _run(("CentralDifferenceLadruno",), dt_bulk, T)
    t_sms, u_sms, sms_div = _run(("CentralDifferenceSMS", DT_TARGET), dt_bulk, T)
    assert cd_div, (
        "plain CentralDifferenceLadruno at the bulk dt=%.4g must DIVERGE (tiny elements "
        "above their stable step); else SMS adds nothing and the test is vacuous" % dt_bulk
    )
    assert not sms_div and math.isfinite(u_sms[-1]), (
        "CentralDifferenceSMS at the bulk dt=%.4g must stay bounded (it raised the tiny "
        "elements' step)" % dt_bulk
    )

    # --- FIDELITY (modal): the ~20%% fictitious mass sits in a LOW-participation zone (near
    #     the clamp), so the fundamental frequency barely moves. This is the honest selling
    #     point: a lot of added mass, almost no structural-frequency cost.
    f1_0 = _f1(after_sms=False)
    f1_s = _f1(after_sms=True)
    df1 = abs(f1_s - f1_0) / f1_0
    assert df1 < 0.01, (
        "fundamental frequency shifted %.3f%% under %.1f%% added mass (expect <1%% for a "
        "low-participation refinement zone)" % (100 * df1, added_pct)
    )

    # --- FIDELITY (response): the SMS run at the bulk dt reproduces the free-end axial
    #     vibration of the fine-dt reference — period and peak amplitude track. (Period/peak
    #     isolate the scaling fidelity from bulk-dt dispersion.)
    t_ref, u_ref, _ = _run(("CentralDifferenceLadruno",), 0.9 * dtcr, T)
    P_ref, P_sms = _period(t_ref, u_ref), _period(t_sms, u_sms)
    assert P_ref and P_sms, "could not extract a vibration period (too few crossings)"
    dP = abs(P_sms - P_ref) / P_ref
    assert dP < 0.03, (
        "free-end period under SMS=%.4f vs fine-dt reference=%.4f (%.2f%%) — structural "
        "dynamics not preserved" % (P_sms, P_ref, 100 * dP)
    )
    pk_ref = max(abs(v) for v in u_ref)
    pk_sms = max(abs(v) for v in u_sms)
    assert 0.85 < pk_sms / pk_ref < 1.15, (
        "free-end peak amplitude ratio SMS/ref=%.3f outside [0.85,1.15] — response not "
        "faithful" % (pk_sms / pk_ref)
    )

    # surface the headline numbers (visible with -s / on failure)
    print("\nT-ZONEB: %d hexes, fine band %d scaled, added mass %.1f%% of model, "
          "df1=%.3f%%, period dP=%.2f%%, peak ratio=%.3f"
          % (len(_HEXES), n_scaled, added_pct, 100 * df1, 100 * dP, pk_sms / pk_ref))
